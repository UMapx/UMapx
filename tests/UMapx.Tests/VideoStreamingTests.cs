using System.Drawing;
using System.Drawing.Imaging;
using System.Net;
using System.Net.Sockets;
using System.Runtime.InteropServices;
using System.Runtime.Versioning;
using System.Text;
using UMapx.Video;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Video")]
[SupportedOSPlatform("windows")]
public class VideoStreamingTests
{
    [Theory]
    [InlineData(false)] [InlineData(true)]
    public async Task MjpegDeliversEveryBufferedFrameInOrderBeforeReadingAgain(bool keepOpen)
    {
        Color[] expected = Enumerable.Range(0, 20).Select(i => (i % 3) switch
        {
            0 => Color.Red, 1 => Color.Blue, _ => Color.Lime
        }).ToArray();
        byte[] payload = Multipart(expected.Select(SmallJpeg));
        var result = await Receive(payload, true, expected.Length, keepOpen: keepOpen);
        Assert.Empty(result.Errors);
        Assert.Equal(expected.Length, result.Count);
        Assert.Equal(expected.Length, result.Frames.Count);
        for (int i = 0; i < expected.Length; i++)
            ImagingAuditTests.Pixel(expected[i], result.Frames[i].Color, 4, false);
    }

    [Fact]
    public async Task MjpegConsumesFramesWithoutSubscribers()
    {
        byte[] payload = Multipart(Enumerable.Repeat(SmallJpeg(Color.Red), 20));
        var result = await Receive(payload, true, 20, subscribe: false);
        Assert.Equal(20, result.Count);
        Assert.Empty(result.Frames);
        // The second request gets a deliberate 503 to terminate this subscriber-free run.
        Assert.Single(result.Errors);
        Assert.Contains("503", result.Errors[0]);
    }

    [Fact]
    public async Task MjpegStopsInsideABatchWhenTheSubscriberRequestsIt()
    {
        byte[] payload = Multipart(Enumerable.Repeat(SmallJpeg(Color.Red), 20));
        var result = await Receive(payload, true, 1);
        Assert.Empty(result.Errors);
        Assert.Equal(1, result.Count);
        Assert.Single(result.Frames);
    }

    [Theory]
    [InlineData(false, false)] [InlineData(false, true)]
    [InlineData(true, false)] [InlineData(true, true)]
    public async Task StreamsDecodeLargeJpegsWithOrWithoutContentLength(bool mjpeg, bool contentLength)
    {
        byte[] jpeg = LargeJpeg();
        Assert.True(jpeg.Length > 3 * 1024 * 1024);
        using var encoded = new MemoryStream(jpeg);
        using var expected = new Bitmap(encoded);
        byte[] payload = mjpeg ? Multipart(new[] { jpeg, SmallJpeg(Color.Blue) }) : jpeg;
        var result = await Receive(payload, mjpeg, mjpeg ? 2 : 1, contentLength: contentLength);
        Assert.Empty(result.Errors);
        Assert.Equal(mjpeg ? 2 : 1, result.Count);
        Assert.Equal(mjpeg ? 2 : 1, result.Frames.Count);
        Assert.Equal(expected.Size, result.Frames[0].Size);
        Assert.Equal(expected.GetPixel(expected.Width / 2, expected.Height / 2), result.Frames[0].Color);
        if (mjpeg) ImagingAuditTests.Pixel(Color.Blue, result.Frames[1].Color, 4, false);
    }

    [Theory]
    [InlineData(0, false)] [InlineData(1, false)] [InlineData(17, false)]
    [InlineData(1024, false)] [InlineData(1048576, false)]
    [InlineData(0, true)] [InlineData(17, true)] [InlineData(1048576, true)]
    public void ParserGrowthPreservesPartialMarkersAndUnreadFrames(int initialCapacity, bool raw)
    {
        byte[] header = { 255, 216, 255 };
        byte[] delimiter = raw ? header : Encoding.ASCII.GetBytes("--frame");
        var data = new List<byte>();
        int[] lengths = { 1023, 1048575, 2097151, 31 };
        foreach (int length in lengths)
        {
            data.AddRange(header);
            data.AddRange(Enumerable.Repeat((byte)42, length - header.Length));
            if (!raw) data.AddRange(delimiter);
        }
        if (raw) data.AddRange(header);
        byte[] payload = data.ToArray();
        var parser = new MJPEGStreamParser(new Boundary(raw ? "" : "--frame") { IsChecked = true }, header, initialCapacity);
        using var input = new ShortReadStream(payload);
        int consumed = 0, frames = 0;
        while (input.Position < input.Length)
        {
            parser.Read(input);
            parser.DetectFrame();
            while (parser.HasFrame)
            {
                int length = lengths[frames++];
                Assert.Equal(header, parser.Content.Take(header.Length));
                Assert.True(parser.Content.AsSpan(header.Length, length - header.Length).IndexOfAnyExcept((byte)42) < 0);
                consumed += length + (raw ? 0 : delimiter.Length);
                parser.RemoveFrame();
                int remaining = (int)input.Position - consumed;
                Assert.Equal(payload.Skip(consumed).Take(remaining), parser.Content.Take(remaining));
                parser.DetectFrame();
            }
        }
        Assert.Equal(4, frames);
        Assert.False(parser.HasFrame);
    }

    private sealed class ShortReadStream(byte[] data) : MemoryStream(data)
    {
        public override int Read(byte[] buffer, int offset, int count) => base.Read(buffer, offset, Math.Min(997, count));
    }

    private sealed record Result(List<(Size Size, Color Color)> Frames, List<string> Errors, int Count);

    private static async Task<Result> Receive(byte[] payload, bool mjpeg, int expectedFrames,
        bool subscribe = true, bool keepOpen = false, bool contentLength = true)
    {
        using var listener = new TcpListener(IPAddress.Loopback, 0);
        listener.Start();
        using var cancellation = new CancellationTokenSource(TimeSpan.FromSeconds(10));
        var release = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
        var server = Task.Run(async () =>
        {
            try
            {
                for (int request = 0; request < 2; request++)
                {
                    using var client = await listener.AcceptTcpClientAsync(cancellation.Token);
                    using var stream = client.GetStream();
                    using var reader = new StreamReader(stream, Encoding.ASCII, leaveOpen: true);
                    while (!string.IsNullOrEmpty(await reader.ReadLineAsync(cancellation.Token))) { }
                    if (request != 0)
                    {
                        await stream.WriteAsync(Encoding.ASCII.GetBytes("HTTP/1.1 503 End of test\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"), cancellation.Token);
                        break;
                    }
                    string type = mjpeg ? "multipart/x-mixed-replace; boundary=frame" : "image/jpeg";
                    string length = contentLength ? $"Content-Length: {payload.Length + (keepOpen ? 1 : 0)}\r\n" : "";
                    await stream.WriteAsync(Encoding.ASCII.GetBytes($"HTTP/1.1 200 OK\r\nContent-Type: {type}\r\n{length}Connection: close\r\n\r\n"), cancellation.Token);
                    await stream.WriteAsync(payload, cancellation.Token);
                    if (keepOpen) await release.Task.WaitAsync(cancellation.Token);
                }
            }
            catch (OperationCanceledException) when (cancellation.IsCancellationRequested) { }
        });
        string url = $"http://127.0.0.1:{((IPEndPoint)listener.LocalEndpoint).Port}/frame";
        using IVideoSource source = mjpeg
            ? new MJPEGStream(url) { RequestTimeout = 2000, Proxy = new WebProxy() }
            : new JPEGStream(url) { RequestTimeout = 2000, Proxy = new WebProxy(), PreventCaching = false };
        var finished = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
        var frames = new List<(Size, Color)>();
        var errors = new List<string>();
        if (subscribe) source.NewFrame += (_, e) =>
        {
            frames.Add((e.Frame.Size, e.Frame.GetPixel(e.Frame.Width / 2, e.Frame.Height / 2)));
            if (frames.Count == expectedFrames) source.SignalToStop();
        };
        source.VideoSourceError += (_, e) => { errors.Add(e.Description); source.SignalToStop(); };
        source.PlayingFinished += (_, _) => finished.TrySetResult();
        try
        {
            source.Start();
            await finished.Task.WaitAsync(TimeSpan.FromSeconds(8));
            source.WaitForStop();
            return new Result(frames, errors, source.FramesReceived);
        }
        finally
        {
            source.SignalToStop();
            release.TrySetResult();
            cancellation.Cancel();
            await server;
            source.WaitForStop();
        }
    }

    private static byte[] Multipart(IEnumerable<byte[]> frames)
    {
        using var output = new MemoryStream();
        foreach (byte[] frame in frames)
        {
            output.Write(Encoding.ASCII.GetBytes("--frame\r\nContent-Type: image/jpeg\r\n\r\n"));
            output.Write(frame);
            output.Write(Encoding.ASCII.GetBytes("\r\n"));
        }
        output.Write(Encoding.ASCII.GetBytes("--frame--\r\n"));
        return output.ToArray();
    }

    private static byte[] SmallJpeg(Color color)
    {
        using var bitmap = new Bitmap(8, 6);
        using (var graphics = Graphics.FromImage(bitmap)) graphics.Clear(color);
        using var output = new MemoryStream();
        bitmap.Save(output, ImageFormat.Jpeg);
        return output.ToArray();
    }

    private static byte[] LargeJpeg()
    {
        using var bitmap = new Bitmap(2048, 2048, PixelFormat.Format24bppRgb);
        var data = bitmap.LockBits(new Rectangle(Point.Empty, bitmap.Size), ImageLockMode.WriteOnly, bitmap.PixelFormat);
        try
        {
            byte[] pixels = new byte[data.Stride * data.Height];
            new Random(123).NextBytes(pixels);
            Marshal.Copy(pixels, 0, data.Scan0, pixels.Length);
        }
        finally { bitmap.UnlockBits(data); }
        using var output = new MemoryStream();
        using var parameters = new EncoderParameters(1);
        parameters.Param[0] = new EncoderParameter(System.Drawing.Imaging.Encoder.Quality, 90L);
        var encoder = ImageCodecInfo.GetImageEncoders().Single(e => e.FormatID == ImageFormat.Jpeg.Guid);
        bitmap.Save(output, encoder, parameters);
        return output.ToArray();
    }
}
