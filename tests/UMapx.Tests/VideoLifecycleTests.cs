using System.Drawing;
using System.Reflection;
using System.Runtime.Versioning;
using UMapx.Video;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Video")]
[SupportedOSPlatform("windows")]
public class VideoLifecycleTests
{
    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void DisposingAnActiveImageSourceStopsItAndPreventsRestart(bool depth)
    {
        using var source = ImageSource(depth);
        source.Start();
        Assert.True(source.IsRunning);
        source.Dispose();
        Assert.False(source.IsRunning);
        Assert.Throws<ObjectDisposedException>(source.Start);
        source.SignalToStop();
        source.Dispose();
    }

    [Theory]
    [InlineData(false, false)] [InlineData(true, false)]
    [InlineData(false, true)] [InlineData(true, true)]
    public void QueuedImageCallbacksDoNotDeliverAfterStopOrIntoANewRun(bool depth, bool restart)
    {
        using var source = ImageSource(depth);
        var oldTimer = StartWithoutAutomaticTicks(source);
        int frames = 0;
        source.NewFrame += (_, _) => frames++;
        source.SignalToStop();
        var currentTimer = restart ? StartWithoutAutomaticTicks(source) : oldTimer;
        Tick(source, oldTimer);
        Assert.Equal(0, frames);
        if (restart)
        {
            Tick(source, currentTimer);
            Assert.Equal(1, frames);
        }
        source.Dispose();
        Tick(source, currentTimer);
        Assert.Equal(restart ? 1 : 0, frames);
    }

    [Theory]
    [InlineData(0)] [InlineData(1)] [InlineData(2)]
    public void ImageFramesAreDisposedWhenFrameOrDepthSubscribersThrow(int mode)
    {
        using var source = ImageSource(mode != 0);
        var timer = StartWithoutAutomaticTicks(source);
        Bitmap? delivered = null;
        var failure = new InvalidOperationException("Subscriber failure");
        source.NewFrame += (_, e) =>
        {
            delivered = e.Frame;
            if (mode != 2) throw failure;
        };
        if (mode == 2) ((VideoImageDepthSource)source).NewDepth += (_, _) => throw failure;
        var exception = Assert.Throws<TargetInvocationException>(() => Tick(source, timer));
        Assert.Same(failure, exception.InnerException);
        Assert.NotNull(delivered);
        Assert.Throws<ArgumentException>(() => delivered.GetPixel(0, 0));
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void DisposingFromAFrameCallbackKeepsTheCurrentFrameAliveButStopsFurtherDepth(bool depth)
    {
        using var source = ImageSource(depth);
        var timer = StartWithoutAutomaticTicks(source);
        Bitmap? delivered = null;
        int depthFrames = 0;
        source.NewFrame += (_, _) => source.Dispose();
        source.NewFrame += (_, e) =>
        {
            delivered = e.Frame;
            Assert.Equal(new Size(2, 2), e.Frame.Size);
        };
        if (depth) ((VideoImageDepthSource)source).NewDepth += (_, _) => depthFrames++;
        Tick(source, timer);
        Assert.False(source.IsRunning);
        Assert.Equal(0, depthFrames);
        Assert.NotNull(delivered);
        Assert.Throws<ArgumentException>(() => delivered.GetPixel(0, 0));
    }

    [Theory]
    [InlineData("jpeg")] [InlineData("mjpeg")] [InlineData("screen")]
    public void DisposedThreadSourcesRejectStart(string kind)
    {
        using IVideoSource source = kind switch
        {
            "jpeg" => new JPEGStream("not a valid URI"),
            "mjpeg" => new MJPEGStream("not a valid URI"),
            _ => new ScreenCaptureStream(new Rectangle(0, 0, 2, 2))
        };
        source.Dispose();
        // Do not start a real screen capture even if the disposed-state guard regresses.
        if (source is ScreenCaptureStream)
            typeof(ScreenCaptureStream).GetField("thread", BindingFlags.Instance | BindingFlags.NonPublic)!
                .SetValue(source, Thread.CurrentThread);
        try { Assert.Throws<ObjectDisposedException>(source.Start); }
        finally
        {
            if (source is ScreenCaptureStream)
                typeof(ScreenCaptureStream).GetField("thread", BindingFlags.Instance | BindingFlags.NonPublic)!
                    .SetValue(source, null);
            else
            {
                source.SignalToStop();
                source.WaitForStop();
            }
        }
        Assert.False(source.IsRunning);
    }

    [Theory]
    [InlineData("jpeg,external")] [InlineData("mjpeg,external")]
    [InlineData("jpeg,callback")] [InlineData("mjpeg,callback")]
    [InlineData("jpeg,concurrent")] [InlineData("mjpeg,concurrent")]
    [InlineData("jpeg,finished")] [InlineData("mjpeg,finished")]
    [InlineData("jpeg,restart")] [InlineData("mjpeg,restart")]
    [InlineData("jpeg,frame-dispose")] [InlineData("mjpeg,frame-dispose")]
    [InlineData("jpeg,frame-throw")] [InlineData("mjpeg,frame-throw")]
    [InlineData("screen,external")]
    public async Task ThreadSourcesStopWithoutCrashingOrDeadlocking(string scenario)
    {
        Assert.Equal("True", await AuditProcess.RunAsync("VideoLifecycle", scenario));
    }

    private static IVideoSource ImageSource(bool depth)
    {
        var resolution = new VideoCapabilities(new Size(2, 2), 1, 1, 32);
        return depth ? new VideoImageDepthSource(new Bitmap(2, 2), new ushort[2, 2])
        {
            VideoResolution = resolution, DepthResolution = resolution
        } : new VideoImageSource(new Bitmap(2, 2)) { VideoResolution = resolution };
    }

    // Exercise queued callbacks deterministically, without racing a real timer tick.
    private static System.Timers.Timer StartWithoutAutomaticTicks(IVideoSource source)
    {
        source.Start();
        var timer = (System.Timers.Timer)source.GetType()
            .GetField("_timer", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(source)!;
        timer.Stop();
        return timer;
    }

    private static void Tick(IVideoSource source, object timer) => source.GetType()
        .GetMethod("OnElapsed", BindingFlags.Instance | BindingFlags.NonPublic)!
        .Invoke(source, new object?[] { timer, null });
}
