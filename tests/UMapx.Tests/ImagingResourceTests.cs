using System.Drawing;
using System.Drawing.Imaging;
using System.Reflection;
using System.Runtime.Versioning;
using UMapx.Core;
using UMapx.Imaging;
using UMapx.Transform;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Imaging")]
[SupportedOSPlatform("windows")]
public class ImagingResourceTests
{
    public static IEnumerable<object[]> TwoImageFilters() => typeof(IBitmapFilter2).Assembly.GetTypes()
        .Where(t => t.Namespace == "UMapx.Imaging" && !t.IsAbstract && typeof(IBitmapFilter2).IsAssignableFrom(t))
        .OrderBy(t => t.Name).Select(t => new object[] { t.Name });

    [Theory, MemberData(nameof(TwoImageFilters))]
    public void BusySecondImageDoesNotLeaveFirstImageLocked(string name)
    {
        using var destination = Bitmap();
        using var source = Bitmap();
        var held = BitmapFormat.Lock32bpp(source);
        try
        {
            Assert.Throws<InvalidOperationException>(() => CreateFilter(name).Apply(destination, source));
            AssertReusable(destination);
            AssertLocked(source);
        }
        finally { BitmapFormat.Unlock(source, held); }
        AssertReusable(source);
    }

    [Theory, MemberData(nameof(TwoImageFilters))]
    public void DisposedSecondImageDoesNotLeaveFirstImageLocked(string name)
    {
        using var destination = Bitmap();
        using var source = Bitmap();
        source.Dispose();
        Assert.Throws<ArgumentException>(() => CreateFilter(name).Apply(destination, source));
        AssertReusable(destination);
    }

    [Theory, MemberData(nameof(TwoImageFilters))]
    public void AliasedImagesDoNotRemainLockedAfterTheSecondLockFails(string name)
    {
        using var image = Bitmap();
        Assert.Throws<InvalidOperationException>(() => CreateFilter(name).Apply(image, image));
        AssertReusable(image);
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void ValidationFailureReleasesBothImages(bool crop)
    {
        using var destination = Bitmap(2);
        using var source = Bitmap(4);
        IBitmapFilter2 filter = crop ? new Crop(0, 0, 3, 3) : new BoxBlur();
        Assert.Throws<ArgumentException>(() => filter.Apply(destination, source));
        AssertReusable(destination);
        AssertReusable(source);
    }

    [Fact]
    public void UserFilterFailurePreservesItsExceptionAndReleasesTheBitmap()
    {
        using var image = Bitmap();
        var failure = new InvalidOperationException("User filter failure");
        var filter = new BitmapFilter(new FailingFilter(failure));
        Assert.Same(failure, Assert.Throws<InvalidOperationException>(() => filter.Apply(image)));
        AssertReusable(image);
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void RebuildFailurePreservesExceptionAndCallerBitmapOwnership(bool borrowedData)
    {
        using var image = Bitmap();
        var failure = new InvalidOperationException("Rebuild failure");
        var filter = new FailingPointAddition(failure);
        if (borrowedData)
        {
            var held = BitmapFormat.Lock32bpp(image);
            try
            {
                Assert.Same(failure, Assert.Throws<InvalidOperationException>(() => filter.Apply(held)));
                AssertLocked(image);
            }
            finally { BitmapFormat.Unlock(image, held); }
        }
        else Assert.Same(failure, Assert.Throws<InvalidOperationException>(() => filter.Apply(image)));
        AssertReusable(image);
    }

    [Theory]
    [InlineData("RGB")] [InlineData("HSB")] [InlineData("HSL")] [InlineData("YCbCr")]
    public void InvalidChannelArraysDoNotLeaveCallerBitmapLocked(string model)
    {
        using var image = Bitmap();
        var channels = new[] { new float[4, 4] };
        var convert = typeof(BitmapMatrix).GetMethod("From" + model, new[] { typeof(float[][,]), typeof(Bitmap) })!;
        var error = Assert.Throws<TargetInvocationException>(() => convert.Invoke(null, new object[] { channels, image }));
        Assert.IsType<IndexOutOfRangeException>(error.InnerException);
        AssertReusable(image);
    }

    [Fact]
    public void StereoDisparityReleasesItsFirstImageIfTheSecondIsBusy()
    {
        using var destination = Bitmap();
        using var source = Bitmap();
        var held = BitmapFormat.Lock32bpp(source);
        try
        {
            Assert.Throws<InvalidOperationException>(() => new StereoDisparity().Apply(destination, source));
            AssertReusable(destination);
            AssertLocked(source);
        }
        finally { BitmapFormat.Unlock(source, held); }
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void MotionDetectorCanRecoverAfterItsStoredFrameWasBusy(bool borrowedData)
    {
        using var image = Bitmap();
        using var detector = new MotionDetector();
        Assert.Equal(0, detector.Apply(image));
        var field = typeof(MotionDetector).GetField("Frame", BindingFlags.Instance | BindingFlags.NonPublic)!;
        var previous = (Bitmap)field.GetValue(detector)!;
        var previousLock = BitmapFormat.Lock32bpp(previous);
        try
        {
            if (borrowedData)
            {
                var held = BitmapFormat.Lock32bpp(image);
                try
                {
                    Assert.Throws<InvalidOperationException>(() => detector.Apply(held));
                    AssertLocked(image);
                }
                finally { BitmapFormat.Unlock(image, held); }
            }
            else Assert.Throws<InvalidOperationException>(() => detector.Apply(image));
            AssertReusable(image);
            Assert.Same(previous, field.GetValue(detector));
            AssertLocked(previous);
        }
        finally { BitmapFormat.Unlock(previous, previousLock); }
        Assert.Equal(0, detector.Apply(image));
    }

    private static IBitmapFilter2 CreateFilter(string name)
    {
        IBitmapFilter2? configured = name switch
        {
            nameof(Crop) => new Crop(0, 0, 4, 4),
            nameof(GaussianBlur) => new GaussianBlur(3, 3),
            nameof(HomomorphicEnhancement) => new HomomorphicEnhancement(3, Space.RGB),
            nameof(KsiContrastEnhancement) => new KsiContrastEnhancement(3, Space.RGB),
            nameof(LocalContrastEnhancement) => new LocalContrastEnhancement(3, Space.RGB),
            nameof(LocalContrastInversion) => new LocalContrastInversion(3, Space.RGB),
            nameof(LocalThreshold) => new LocalThreshold(3, Space.RGB),
            nameof(MaskColorFilter) => new MaskColorFilter(Color.Red),
            nameof(Operation) => new Operation(1, 1),
            nameof(Rotate) => new Rotate(0),
            nameof(ShadowsHighlightsCorrection) => new ShadowsHighlightsCorrection(3, Space.RGB),
            nameof(SingleScaleRetinex) => new SingleScaleRetinex(3, Space.RGB),
            nameof(StereoAnaglyph) => new StereoAnaglyph(AnaglyphMode.Color),
            _ => null
        };
        if (configured != null) return configured;
        var type = typeof(IBitmapFilter2).Assembly.GetType("UMapx.Imaging." + name)!;
        var constructor = type.GetConstructors().OrderBy(c => c.GetParameters().Length)
            .First(c => c.GetParameters().All(p => p.HasDefaultValue));
        return (IBitmapFilter2)constructor.Invoke(constructor.GetParameters().Select(p => p.DefaultValue).ToArray());
    }

    private static Bitmap Bitmap(int size = 4) => new(size, size, PixelFormat.Format32bppArgb);
    private static void AssertReusable(Bitmap image)
    {
        var data = BitmapFormat.Lock32bpp(image);
        try { Assert.Equal(image.Width, data.Width); }
        finally { BitmapFormat.Unlock(image, data); }
    }
    private static void AssertLocked(Bitmap image) => Assert.Throws<InvalidOperationException>(() => AssertReusable(image));

    private sealed class FailingPointAddition : PointAddition
    {
        private readonly Exception failure;
        internal FailingPointAddition(Exception failure)
        {
            this.failure = failure;
            rebuild = true;
        }
        protected override void Rebuild() => throw failure;
    }
    private sealed class FailingFilter(Exception failure) : IFilter
    {
        public void Apply(float[] data) => throw failure;
        public void Apply(float[,] data) => throw failure;
        public void Apply(Complex32[] data) => throw failure;
        public void Apply(Complex32[,] data) => throw failure;
    }
}
