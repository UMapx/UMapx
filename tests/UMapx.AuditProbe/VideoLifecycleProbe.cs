using System.Drawing;
using System.Drawing.Imaging;
using System.Net;
using System.Net.Sockets;
using System.Reflection;
using System.Runtime.Versioning;
using System.Text;
using UMapx.Video;

namespace UMapx.AuditProbe;

internal static class VideoLifecycleProbe
{
    internal static bool Run(string argument)
    {
        string[] parts = argument.Split(',');
        if (parts[1].StartsWith("blocked-", StringComparison.Ordinal))
            return BlockedNetwork(parts[0] == "mjpeg", parts[1]);
        if (parts[1].StartsWith("frame-", StringComparison.Ordinal))
        {
            if (!OperatingSystem.IsWindows()) throw new PlatformNotSupportedException();
            return FrameDisposal(parts[0] == "mjpeg", parts[1] == "frame-throw");
        }
        if (parts[0] == "screen")
        {
            if (!OperatingSystem.IsWindows()) throw new PlatformNotSupportedException();
            return ScreenDisposal();
        }
        string mode = parts[1];
        using IVideoSource source = parts[0] == "jpeg"
            ? new JPEGStream("not a valid URI") : new MJPEGStream("not a valid URI");
        using var entered = new ManualResetEventSlim();
        using var resume = new ManualResetEventSlim();
        using var finished = new ManualResetEventSlim();
        Thread? worker = null;
        int errors = 0, finishes = 0;
        source.VideoSourceError += (_, _) =>
        {
            worker = Thread.CurrentThread;
            Interlocked.Increment(ref errors);
            entered.Set();
            switch (mode)
            {
                case "callback": source.Dispose(); break;
                case "finished":
                case "restart": source.SignalToStop(); break;
                default:
                    Require(resume.Wait(TimeSpan.FromSeconds(2)), "Caller did not resume the callback");
                    if (mode == "concurrent") source.Dispose();
                    break;
            }
        };
        source.PlayingFinished += (_, reason) =>
        {
            Require(reason == ReasonToFinishPlaying.StoppedByUser, "Unexpected stop reason");
            if (mode == "finished") source.Dispose();
            Interlocked.Increment(ref finishes);
            finished.Set();
        };
        for (int run = 0; run < (mode == "restart" ? 2 : 1); run++)
        {
            entered.Reset();
            finished.Reset();
            source.Start();
            Require(entered.Wait(TimeSpan.FromSeconds(2)), "Worker did not start");
            if (mode is "external" or "concurrent")
            {
                var stop = StopEvent(source);
                var disposal = Task.Run(source.Dispose);
                try
                {
                    Require(stop.WaitOne(TimeSpan.FromSeconds(2)), "Dispose did not request a stop");
                    Require(!disposal.IsCompleted, "Dispose returned while a callback was still running");
                }
                finally { resume.Set(); }
                Require(disposal.Wait(TimeSpan.FromSeconds(2)), "Dispose deadlocked");
            }
            Require(finished.Wait(TimeSpan.FromSeconds(2)), "Worker did not notify completion");
            Require(worker!.Join(TimeSpan.FromSeconds(2)), "Worker did not terminate");
            source.WaitForStop();
            Require(!source.IsRunning, "Source still reports running");
        }
        source.Dispose();
        source.Dispose();
        source.SignalToStop();
        Require(errors == (mode == "restart" ? 2 : 1), "Unexpected error callback count");
        Require(finishes == errors, "Missing or duplicate completion notification");
        try { source.Start(); }
        catch (ObjectDisposedException) { return true; }
        throw new Exception("Disposed source restarted");
    }

    [SupportedOSPlatform("windows")]
    private static bool ScreenDisposal()
    {
        using var source = new ScreenCaptureStream(new Rectangle(0, 0, 2, 2));
        using var entered = new ManualResetEventSlim();
        using var stop = new ManualResetEvent(false);
        var type = source.GetType();
        type.GetField("stopEvent", BindingFlags.Instance | BindingFlags.NonPublic)!.SetValue(source, stop);
        // Enter the real worker only after stop is signalled, so no screen is captured.
        var worker = new Thread(() =>
        {
            entered.Set();
            Require(stop.WaitOne(TimeSpan.FromSeconds(2)), "Dispose did not stop screen capture");
            type.GetMethod("WorkerThread", BindingFlags.Instance | BindingFlags.NonPublic)!.Invoke(source, null);
        });
        type.GetField("thread", BindingFlags.Instance | BindingFlags.NonPublic)!.SetValue(source, worker);
        worker.Start();
        Require(entered.Wait(TimeSpan.FromSeconds(2)), "Worker did not start");
        source.Dispose();
        Require(worker.Join(TimeSpan.FromSeconds(2)), "Screen worker did not terminate");
        return !source.IsRunning;
    }

    [SupportedOSPlatform("windows")]
    private static bool FrameDisposal(bool mjpeg, bool throws)
    {
        using var input = new Bitmap(2, 2);
        using var encoded = new MemoryStream();
        input.Save(encoded, ImageFormat.Jpeg);
        byte[] payload = encoded.ToArray();
        if (mjpeg)
            payload = Encoding.ASCII.GetBytes("--frame\r\nContent-Type: image/jpeg\r\n\r\n")
                .Concat(payload).Concat(Encoding.ASCII.GetBytes("\r\n--frame--\r\n")).ToArray();
        string contentType = mjpeg ? "multipart/x-mixed-replace; boundary=frame" : "image/jpeg";
        var listener = new TcpListener(IPAddress.Loopback, 0);
        listener.Start();
        try
        {
            var server = Task.Run(async () =>
            {
                using var client = await listener.AcceptTcpClientAsync();
                using var stream = client.GetStream();
                using var reader = new StreamReader(stream, Encoding.ASCII, leaveOpen: true);
                while (!string.IsNullOrEmpty(await reader.ReadLineAsync())) { }
                byte[] header = Encoding.ASCII.GetBytes($"HTTP/1.1 200 OK\r\nContent-Type: {contentType}\r\nContent-Length: {payload.Length}\r\nConnection: close\r\n\r\n");
                await stream.WriteAsync(header);
                await stream.WriteAsync(payload);
            });
            string url = $"http://127.0.0.1:{((IPEndPoint)listener.LocalEndpoint).Port}/frame";
            using IVideoSource source = mjpeg
                ? new MJPEGStream(url) { RequestTimeout = 1000, Proxy = new WebProxy() }
                : new JPEGStream(url) { RequestTimeout = 1000, Proxy = new WebProxy() };
            using var finished = new ManualResetEventSlim();
            Bitmap? delivered = null;
            Thread? worker = null;
            int errors = 0;
            string? error = null;
            source.NewFrame += (_, e) =>
            {
                worker = Thread.CurrentThread;
                delivered = e.Frame;
                if (throws) throw new InvalidOperationException("Subscriber failure");
                source.Dispose();
                Require(e.Frame.Size == new Size(2, 2), "Current frame disposed inside its callback");
            };
            source.VideoSourceError += (_, e) =>
            {
                error = e.Description;
                Interlocked.Increment(ref errors);
                source.Dispose();
            };
            source.PlayingFinished += (_, _) => finished.Set();
            source.Start();
            Require(finished.Wait(TimeSpan.FromSeconds(3)), "Frame source did not finish");
            Require(worker != null && worker.Join(TimeSpan.FromSeconds(1)), "No frame delivered or worker still running: " + error);
            Require(server.Wait(TimeSpan.FromSeconds(1)), "Local server did not finish");
            Require(errors == (throws ? 1 : 0), "Unexpected source error");
            Require(!source.IsRunning, "Frame source still reports running");
            try { delivered!.GetPixel(0, 0); }
            catch (ArgumentException) { return true; }
            finally { delivered?.Dispose(); }
            throw new Exception("Delivered bitmap was not disposed");
        }
        finally { listener.Stop(); }
    }

    private static bool BlockedNetwork(bool mjpeg, string scenario)
    {
        string[] parts = scenario.Split('-');
        string mode = parts[2];
        using var listener = new TcpListener(IPAddress.Loopback, 0);
        listener.Start();
        string url = $"http://127.0.0.1:{((IPEndPoint)listener.LocalEndpoint).Port}/frame";
        using IVideoSource source = mjpeg
            ? new MJPEGStream(url) { RequestTimeout = Timeout.Infinite, Proxy = new WebProxy() }
            : new JPEGStream(url) { RequestTimeout = Timeout.Infinite, Proxy = new WebProxy(), PreventCaching = false };
        int errors = 0, frames = 0, finishes = 0;
        ReasonToFinishPlaying? reason = null;
        source.NewFrame += (_, _) => Interlocked.Increment(ref frames);
        source.VideoSourceError += (_, _) => Interlocked.Increment(ref errors);
        source.PlayingFinished += (_, value) => { reason = value; Interlocked.Increment(ref finishes); };
        int runs = mode == "restart" ? 2 : 1;
        for (int run = 0; run < runs; run++)
        {
            using var entered = new ManualResetEventSlim();
            using var cancellation = new CancellationTokenSource(TimeSpan.FromSeconds(3));
            var release = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
            var server = Task.Run(async () =>
            {
                try
                {
                    using var client = await listener.AcceptTcpClientAsync(cancellation.Token);
                    using var stream = client.GetStream();
                    using var reader = new StreamReader(stream, Encoding.ASCII, leaveOpen: true);
                    while (!string.IsNullOrEmpty(await reader.ReadLineAsync(cancellation.Token))) { }
                    string type = mjpeg ? "multipart/x-mixed-replace; boundary=frame" : "image/jpeg";
                    await stream.WriteAsync(Encoding.ASCII.GetBytes($"HTTP/1.1 200 OK\r\nContent-Type: {type}\r\nContent-Length: 1000\r\nConnection: close\r\n\r\n"), cancellation.Token);
                    // One byte proves the worker entered the body-reading loop; no complete frame follows.
                    await stream.WriteAsync(new byte[] { 255 }, cancellation.Token);
                    entered.Set();
                    await release.Task.WaitAsync(cancellation.Token);
                }
                catch (OperationCanceledException) when (cancellation.IsCancellationRequested) { }
            });
            Task? stopping = null;
            try
            {
                source.Start();
                var worker = (Thread)source.GetType().GetField(mjpeg ? "_thread" : "thread", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(source)!;
                Require(entered.Wait(TimeSpan.FromSeconds(2)), "Server did not receive the request");
                Require(SpinWait.SpinUntil(() => source.BytesReceived > 0, 1000), "Worker did not read the partial body");
                Require(SpinWait.SpinUntil(() => (worker.ThreadState & ThreadState.WaitSleepJoin) != 0, 1000), "Worker did not enter a blocking network operation");
                Action stopAndWait = () => { source.SignalToStop(); source.WaitForStop(); };
                stopping = mode == "concurrent"
                    ? Task.WhenAll(Task.Run(source.Dispose), Task.Run(source.Dispose), Task.Run(stopAndWait))
                    : Task.Run(mode == "dispose" ? source.Dispose : stopAndWait);
                Require(stopping.Wait(TimeSpan.FromSeconds(1)), "Stopping waited for the stalled server");
                Require(!release.Task.IsCompleted, "Server was released before the stop completed");
                Require(!source.IsRunning && !worker.IsAlive, "Worker remained alive after stopping");
                Require(finishes == run + 1 && reason == ReasonToFinishPlaying.StoppedByUser, "Missing or duplicate stop notification");
                Require(errors == 0, "Cancellation was reported as a source error");
                Require(frames == 0, "An incomplete frame was delivered");
            }
            finally
            {
                release.TrySetResult();
                cancellation.Cancel();
                Require(server.Wait(TimeSpan.FromSeconds(2)), "Server did not finish");
                source.SignalToStop();
                source.WaitForStop();
                if (stopping != null) Require(stopping.Wait(TimeSpan.FromSeconds(1)), "Stop did not finish after server cleanup");
            }
        }
        source.Dispose();
        source.Dispose();
        return true;
    }

    private static ManualResetEvent StopEvent(IVideoSource source) => (ManualResetEvent)source.GetType()
        .GetField(source is MJPEGStream ? "_stopEvent" : "stopEvent", BindingFlags.Instance | BindingFlags.NonPublic)!
        .GetValue(source)!;

    private static void Require(bool condition, string message)
    {
        if (!condition) throw new Exception(message);
    }
}
