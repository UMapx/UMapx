using System.Drawing;
using System.Drawing.Imaging;
using System.Globalization;
using System.Net;
using System.Runtime.Versioning;
using System.Text;
using UMapx.Video;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Video")]
[SupportedOSPlatform("windows")]
public class VideoRepairTests
{
    private sealed class Response(string contentType) : WebResponse
    {
        public override string ContentType { get; set; } = contentType;
    }

    [Theory]
    [InlineData("MULTIPART/X-MIXED-REPLACE; BOUNDARY=AaB03x", "AaB03x")]
    [InlineData(" multipart/mixed ; boundary = frame ", "frame")]
    [InlineData("multipart/x-mixed-replace;\tBoundary\t=\tframe", "frame")]
    [InlineData("multipart/x-mixed-replace; boundary=\"part:17\"", "part:17")]
    [InlineData("multipart/x-mixed-replace; boundary=\"Aa B\"; charset=utf-8", "Aa B")]
    [InlineData("multipart/x-mixed-replace; charset=utf-8; boundary=frame", "frame")]
    [InlineData("multipart/x-mixed-replace; boundary=frame; name=video", "frame")]
    [InlineData("multipart/x-mixed-replace; note=\"x; boundary=wrong\"; boundary=Right", "Right")]
    [InlineData("multipart/x-mixed-replace; x-boundary=wrong; boundary=Right", "Right")]
    [InlineData("multipart/x-mixed-replace; x-boundary=wrong", "")]
    [InlineData("multipart/x-mixed-replace", "")]
    [InlineData("multipart/x-mixed-replace; note=\"escaped \\\"quote; boundary=wrong\"; boundary=Right", "Right")]
    [InlineData("multipart/x-mixed-replace;\r\n boundary=frame", "frame")]
    [InlineData("multipart/x-mixed-replace;\r\n\tboundary=frame;\r\n name=video", "frame")]
    [InlineData("multipart/x-mixed-replace; boundary=\"part\\:17\"", "part:17")]
    [InlineData("multipart/x-mixed-replace; boundary=\"\"", "")]
    [InlineData("APPLICATION/OCTET-STREAM; charset=binary", "")]
    [InlineData("application/octet-stream; boundary=ignored", "")]
    [InlineData("Multipart/Mixed; Boundary=\"--frame\"", "--frame")]
    [InlineData("multipart/x-mixed-replace; charset=utf-8; Boundary=\"--Aa:42\"; name=video", "--Aa:42")]
    public void BoundaryParsingSeparatesMimeSyntaxFromTheCaseSensitiveValue(string contentType, string expected)
    {
        using var response = new Response(contentType);
        var boundary = Boundary.FromResponse(response);
        Assert.Equal(expected, boundary.Content);
        Assert.Equal(Encoding.ASCII.GetBytes(expected), (byte[])boundary);
        Assert.Equal(expected.Length != 0, boundary.HasValue);
        Assert.False(boundary.IsChecked);
        Assert.Equal(expected.Length == 0, boundary.IsValid);
    }

    [Theory]
    [InlineData("text/html")]
    [InlineData("image/jpeg")]
    [InlineData("multipart/related; boundary=frame")]
    [InlineData("multipart/mixed-up; boundary=frame")]
    [InlineData("application/octet-stream-extra")]
    [InlineData("multipartish/mixed; boundary=frame")]
    [InlineData("multipart/x-mixed-replace; boundary=\"unterminated")]
    [InlineData("multipart/x-mixed-replace; boundary")]
    [InlineData("multipart/x-mixed-replace;\r\nboundary=frame")]
    [InlineData("multipart/x-mixed-replace;\n boundary=frame")]
    [InlineData("multipart/x-mixed-replace; boundary=first; BOUNDARY=second")]
    [InlineData("multipart/x-mixed-replace; boundary=\"\"; boundary=frame")]
    [InlineData("")]
    [InlineData(" ")]
    public void InvalidOrUnsupportedContentTypesAreRejected(string contentType)
    {
        using var response = new Response(contentType);
        Assert.Throws<ArgumentException>(() => Boundary.FromResponse(response));
    }

    [Fact]
    public void MissingResponseIsRejected()
    {
        Assert.Throws<ArgumentNullException>(() => Boundary.FromResponse(null!));
        using var response = new Response(null!);
        Assert.Throws<ArgumentException>(() => Boundary.FromResponse(response));
    }

    [Theory]
    [InlineData("tr-TR")]
    [InlineData("en-US")]
    public void MimeNameComparisonIsIndependentOfCurrentCulture(string culture)
    {
        var original = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo(culture);
            using var response = new Response("MULTIPART/X-MIXED-REPLACE; BOUNDARY=MiXeD");
            Assert.Equal("MiXeD", Boundary.FromResponse(response).Content);
        }
        finally { CultureInfo.CurrentCulture = original; }
    }

    public static IEnumerable<object[]> FrameCases()
    {
        foreach (int mode in new[] { 0, 1, 2 })
            foreach (int chunk in new[] { 1, 2, 3, 4, 5, 7, 17, 63, 1024 })
                yield return new object[] { mode, chunk };
    }

    [Theory]
    [MemberData(nameof(FrameCases))]
    public void ActualJpegFramesSurvivePartialReadsAndBoundaryCorrection(int mode, int chunk)
    {
        byte[] header = { 255, 216, 255 };
        Color[] expected = { Color.Red, Color.Blue, Color.Lime };
        var data = new List<byte>();
        foreach (var color in expected)
        {
            if (mode != 0) data.AddRange(Encoding.ASCII.GetBytes("--frame\r\nContent-Type: image/jpeg\r\n\r\n"));
            data.AddRange(Jpeg(color));
            if (mode != 0) data.AddRange(Encoding.ASCII.GetBytes("\r\n"));
        }
        // Raw streams delimit a frame with the next JPEG header, including the last tested frame.
        data.AddRange(mode == 0 ? header : Encoding.ASCII.GetBytes("--frame--\r\n"));
        using var response = new Response(mode == 0 ? "Application/Octet-Stream" :
            "Multipart/X-Mixed-Replace; charset=binary; Boundary=" + (mode == 1 ? "frame" : "--frame"));
        var boundary = Boundary.FromResponse(response);
        var parser = new MJPEGStreamParser(boundary, header, 4096);
        using var input = new PartitionStream(data.ToArray(), Enumerable.Repeat(chunk, data.Count).ToArray());
        var actual = new List<Color>();
        parser.DetectFrame();
        Assert.False(parser.HasFrame);
        while (input.Position < input.Length)
        {
            parser.Read(input);
            if (boundary.HasValue && !boundary.IsChecked) boundary.FixMalformedBoundary(parser);
            if (!boundary.IsValid) continue;
            parser.DetectFrame();
            while (parser.HasFrame)
            {
                Assert.True(actual.Count < expected.Length, "A frame was duplicated.");
                using var frame = parser.GetFrame();
                Assert.Equal(new Size(8, 6), frame.Size);
                actual.Add(frame.GetPixel(4, 3));
                parser.RemoveFrame();
                parser.DetectFrame();
            }
        }
        Assert.Equal(expected.Length, actual.Count);
        for (int i = 0; i < expected.Length; i++) ImagingAuditTests.Pixel(expected[i], actual[i], 4, false);
        if (mode != 0) Assert.Equal("--frame", boundary.Content);
        Assert.False(parser.HasFrame);
    }

    [Theory]
    [InlineData(0)] [InlineData(1)] [InlineData(2)] [InlineData(7)] [InlineData(70)]
    public void MarkerSearchRetainsEveryPossiblePartialHeaderAndDelimiter(int boundaryLength)
    {
        string boundary = new('~', boundaryLength);
        byte[] data = MarkerStream(boundary, out int[] removals);
        // Every possible split includes splits inside both headers and delimiters.
        for (int split = 1; split < data.Length; split++)
            AssertMarkerFrames(data, boundary, removals, new[] { split, data.Length - split });
        for (int chunk = 1; chunk <= 9; chunk++)
            AssertMarkerFrames(data, boundary, removals, Enumerable.Repeat(chunk, data.Length).ToArray());
    }

    /// <summary>Checks detection and exact unread bytes against known frame offsets.</summary>
    /// <param name="data">Synthetic frame bytes, independent of JPEG decoding.</param>
    /// <param name="boundary">Delimiter, or empty for JPEG-header delimiters.</param>
    /// <param name="removals">Absolute offsets after removing each complete frame.</param>
    /// <param name="chunks">Transport partition lengths.</param>
    private static void AssertMarkerFrames(byte[] data, string boundary, int[] removals, int[] chunks)
    {
        var parser = new MJPEGStreamParser(new Boundary(boundary) { IsChecked = true }, new byte[] { 255, 216, 255 }, 4096);
        using var input = new PartitionStream(data, chunks);
        int frames = 0;
        parser.DetectFrame();
        Assert.False(parser.HasFrame);
        while (input.Position < input.Length)
        {
            parser.Read(input);
            parser.DetectFrame();
            bool detected = parser.HasFrame;
            parser.DetectFrame();
            Assert.Equal(detected, parser.HasFrame);
            while (parser.HasFrame)
            {
                Assert.True(frames < removals.Length, "A marker was detected twice.");
                int offset = removals[frames++];
                int available = (int)input.Position - offset;
                Assert.True(available >= 0, "A frame was detected before its complete delimiter arrived.");
                parser.RemoveFrame();
                Assert.Equal(data.Skip(offset).Take(available), parser.Content.Take(available));
                parser.DetectFrame();
            }
        }
        Assert.Equal(removals.Length, frames);
    }

    /// <summary>Creates two frames with misleading partial headers and known removal offsets.</summary>
    /// <param name="boundary">Delimiter string, or empty to use the next header.</param>
    /// <param name="removals">Offsets of unread bytes after each removal.</param>
    /// <returns>Stream bytes ending with a delimiter for the second frame.</returns>
    private static byte[] MarkerStream(string boundary, out int[] removals)
    {
        byte[] header = { 255, 216, 255 };
        var data = new List<byte> { 1, 255, 216, 2, 255 };
        removals = new int[2];
        for (int i = 0; i < 2; i++)
        {
            if (boundary.Length == 0 && i > 0) removals[i - 1] = data.Count;
            data.AddRange(header);
            data.AddRange(new byte[] { (byte)(11 + i), 12, 255, 216, 13, 255, 14 });
            if (boundary.Length != 0)
            {
                data.AddRange(Encoding.ASCII.GetBytes(boundary));
                removals[i] = data.Count;
            }
        }
        if (boundary.Length == 0)
        {
            removals[1] = data.Count;
            data.AddRange(header);
        }
        return data.ToArray();
    }

    /// <summary>Encodes a solid-color JPEG using the platform encoder.</summary>
    /// <param name="color">Expected decoded color, subject to JPEG quantization.</param>
    /// <returns>Encoded image bytes.</returns>
    private static byte[] Jpeg(Color color)
    {
        using var bitmap = new Bitmap(8, 6);
        using (var graphics = Graphics.FromImage(bitmap)) graphics.Clear(color);
        using var stream = new MemoryStream();
        bitmap.Save(stream, ImageFormat.Jpeg);
        return stream.ToArray();
    }

    /// <summary>Returns prescribed transport partitions while honoring the caller's read limit.</summary>
    private sealed class PartitionStream(byte[] data, int[] chunks) : MemoryStream(data)
    {
        private int index;
        private int remaining = chunks[0];

        public override int Read(byte[] buffer, int offset, int count)
        {
            int read = base.Read(buffer, offset, Math.Min(count, remaining));
            remaining -= read;
            if (remaining == 0) remaining = ++index < chunks.Length ? chunks[index] : int.MaxValue;
            return read;
        }
    }
}
