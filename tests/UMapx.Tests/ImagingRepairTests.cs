using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
using System.Runtime.Versioning;
using UMapx.Core;
using UMapx.Imaging;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Imaging")]
[SupportedOSPlatform("windows")]
public class ImagingRepairTests
{
    public static IEnumerable<object[]> TableCases()
    {
        foreach (string name in new[] { "Gamma", "Shift", "Bin", "Exposure", "Sin", "Cos", "Log", "Add", "Contrast", "LogContrast", "Invert", "Equalize", "Linear", "Levels" })
            foreach (int length in new[] { 0, 1, 2, 3, 17, 256 })
                yield return new object[] { name, length };
    }

    [Theory]
    [MemberData(nameof(TableCases))]
    public void CorrectionTablesSampleIndependentFunctionsIncludingBothEndpoints(string name, int length)
    {
        float[] actual = name switch
        {
            "Gamma" => Intensity.Gamma(2, length), "Shift" => Intensity.Shift(.25f, length),
            "Bin" => Intensity.Bin(.5f, length), "Exposure" => Intensity.Exposure(127.5f, length),
            "Sin" => Intensity.Sin(.125f, length), "Cos" => Intensity.Cos(.125f, length),
            "Log" => Intensity.Log(2, 0, length), "Add" => Intensity.Add(.125f, length),
            "Contrast" => Intensity.Contrast(.5f, length), "LogContrast" => Intensity.LogContrast(2, length),
            "Invert" => Intensity.Invert(length), "Equalize" => Intensity.Equalize(.25f, .75f, length),
            "Linear" => Intensity.Linear(.75f, .25f, .125f, length),
            _ => Intensity.Levels(.25f, .75f, .125f, .875f, length)
        };
        Assert.Equal(length, actual.Length);
        for (int i = 0; i < length; i++)
        {
            double x = length == 1 ? 0 : i / (double)(length - 1);
            double expected = name switch
            {
                "Gamma" => x * x, "Shift" => Math.Sqrt(x), "Bin" => x > .5 ? 1 : 0,
                "Exposure" => 1 - Math.Exp(-2 * x), "Sin" => .5 * Math.Sin(Math.PI * (x - .5)) + .625,
                "Cos" => .5 * Math.Cos(Math.PI * (x - 1)) + .625, "Log" => Math.Log2(1 + 2 * x),
                "Add" => x + .125, "Contrast" => 1.5 * x - .25,
                "LogContrast" => x <= .5 ? 2 * x * x : 1 - 2 * (1 - x) * (1 - x),
                "Invert" => 1 - x, "Equalize" => Math.Clamp((x - .25) * 2, 0, 1),
                "Linear" => (x - .25) * 2 + .125,
                _ => x <= .25 ? .125 : x >= .75 ? .875 : .125 + (x - .25) * 1.5
            };
            Close(expected, actual[i], 3e-7, 3e-6);
        }
    }

    [Theory]
    [InlineData("Brightness")]
    [InlineData("Contrast")]
    [InlineData("Gamma")]
    [InlineData("Linear")]
    [InlineData("Levels")]
    [InlineData("Shift")]
    [InlineData("Transparency")]
    public void NeutralCorrectionsPreserveEveryByteAcrossRepeatedApplications(string name)
    {
        using var expected = CreateImage(256, 1, (x, y) => Color.FromArgb(x, x, 255 - x, (x * 73) % 256));
        using var actual = (Bitmap)expected.Clone();
        IBitmapFilter filter = name switch
        {
            "Brightness" => new BrightnessCorrection(0, Space.RGB), "Contrast" => new ContrastCorrection(0, Space.RGB),
            "Gamma" => new GammaCorrection(1, Space.RGB), "Linear" => new LinearCorrection(0, Space.RGB),
            "Levels" => new LevelsCorrection(new RangeFloat(0, 1), new RangeFloat(0, 1), Space.RGB),
            "Shift" => new ShiftCorrection(0, Space.RGB), _ => new TransparencyCorrection(0)
        };
        for (int repeat = 0; repeat < 5; repeat++) filter.Apply(actual);
        ImagingAuditTests.Same(expected, actual);
    }

    [Fact]
    public void InversionIsAnExactByteComplementAndAnInvolution()
    {
        using var original = CreateImage(256, 1, (x, y) => Color.FromArgb(x, x, 255 - x, (x * 73) % 256));
        using var actual = (Bitmap)original.Clone();
        var filter = new InvertChannels(Space.RGB);
        filter.Apply(actual);
        for (int x = 0; x < 256; x++)
        {
            var p = original.GetPixel(x, 0);
            ImagingAuditTests.Pixel(Color.FromArgb(p.A, 255 - p.R, 255 - p.G, 255 - p.B), actual.GetPixel(x, 0));
        }
        filter.Apply(actual);
        ImagingAuditTests.Same(original, actual);
    }

    public static IEnumerable<object[]> PixelEquationCases()
    {
        foreach (string name in new[] { "Gamma", "Brightness", "Contrast", "Linear", "Shift", "Transparency" })
            foreach (float value in name == "Shift" ? new[] { -.25f, .125f, .25f } : new[] { .25f, .5f, .75f })
                yield return new object[] { name, value };
    }

    [Theory]
    [MemberData(nameof(PixelEquationCases))]
    public void CorrectionsMatchIndependentByteEquations(string name, float value)
    {
        using var original = CreateImage(256, 1, (x, y) => Color.FromArgb(x, x, 255 - x, (x * 73) % 256));
        using var actual = (Bitmap)original.Clone();
        IBitmapFilter filter = name switch
        {
            "Gamma" => new GammaCorrection(value, Space.RGB), "Brightness" => new BrightnessCorrection(value, Space.RGB),
            "Contrast" => new ContrastCorrection(value, Space.RGB), "Linear" => new LinearCorrection(value, Space.RGB),
            "Shift" => new ShiftCorrection(value, Space.RGB), _ => new TransparencyCorrection(value)
        };
        filter.Apply(actual);
        int Map(int v)
        {
            double x = v / 255.0;
            double result = name switch
            {
                "Gamma" => Math.Pow(x, value), "Brightness" or "Linear" => x + value / 2.0,
                "Contrast" => (x - .5) * (1 + value) + .5,
                "Shift" => Math.Pow(x, Math.Log(.5) / Math.Log(.5 - value)),
                _ => x + value
            };
            return QuantizeByte(255 * result);
        }
        for (int x = 0; x < 256; x++)
        {
            var p = original.GetPixel(x, 0);
            var expected = name == "Transparency" ? Color.FromArgb(Map(p.A), p.R, p.G, p.B)
                : Color.FromArgb(p.A, Map(p.R), Map(p.G), Map(p.B));
            ImagingAuditTests.Pixel(expected, actual.GetPixel(x, 0), 1);
        }
    }

    public static IEnumerable<object[]> StrideCases()
    {
        foreach (string filter in new[] { "RGB", "Grayscale", "Transparency", "Diffusion" })
            foreach (var size in new[] { (1, 1), (1, 7), (7, 1), (3, 5), (7, 4) })
                foreach (int padding in new[] { 0, 4, 20 }) foreach (bool negative in new[] { false, true })
                    yield return new object[] { filter, size.Item1, size.Item2, padding, negative };
    }

    [Theory]
    [MemberData(nameof(StrideCases))]
    public void SignedStrideMatchesScalarPixelsAndProtectsEveryNonpixelByte(string name, int width, int height, int padding, bool negative)
    {
        int pitch = width * 4 + padding, step = negative ? -pitch : pitch;
        var buffer = Enumerable.Repeat((byte)173, 64 + 3 * height * pitch).ToArray();
        int origin = 32 + height * pitch + (negative ? (height - 1) * pitch : 0);
        var active = new bool[buffer.Length];
        var expected = new int[height, width, 4];
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++)
        {
            int k = origin + y * step + 4 * x;
            int[] pixel = { (180 + 23 * x + 17 * y) % 256, (90 + 11 * x + 41 * y) % 256, (30 + 31 * x + 7 * y) % 256, (13 + x * 47 + y * 19) % 256 };
            for (int c = 0; c < 4; c++) { buffer[k + c] = (byte)pixel[c]; active[k + c] = true; expected[y, x, c] = pixel[c]; }
        }
        if (name == "Diffusion") ReferenceFloydSteinberg(expected, width, height);
        else for (int y = 0; y < height; y++) for (int x = 0; x < width; x++)
        {
            if (name == "RGB")
            {
                expected[y, x, 0] = QuantizeByte(expected[y, x, 0] - 15);
                expected[y, x, 1] = QuantizeByte(expected[y, x, 1] + 23);
                expected[y, x, 2] = QuantizeByte(expected[y, x, 2] + 12);
            }
            else if (name == "Transparency") expected[y, x, 3] = QuantizeByte(expected[y, x, 3] + .25 * 255);
            else
            {
                int gray = QuantizeByte(.25 * expected[y, x, 0] + .5 * expected[y, x, 1] + .25 * expected[y, x, 2]);
                for (int c = 0; c < 3; c++) expected[y, x, c] = gray;
            }
        }
        IBitmapFilter filter = name switch
        {
            "RGB" => new RGBFilter(12, 23, -15), "Grayscale" => new Grayscale(.25f, .5f, .25f),
            "Transparency" => new TransparencyCorrection(.25f), _ => ErrorDiffusionDithering.FloydSteinberg
        };
        var handle = GCHandle.Alloc(buffer, GCHandleType.Pinned);
        try
        {
            filter.Apply(new BitmapData { Width = width, Height = height, Stride = step,
                PixelFormat = PixelFormat.Format32bppArgb, Scan0 = IntPtr.Add(handle.AddrOfPinnedObject(), origin) });
        }
        finally { handle.Free(); }
        for (int i = 0; i < buffer.Length; i++) if (!active[i]) Assert.True(buffer[i] == 173, $"Nonpixel byte {i} was modified.");
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++) for (int c = 0; c < 4; c++)
            Assert.Equal(expected[y, x, c], buffer[origin + y * step + 4 * x + c]);
    }

    public static IEnumerable<object[]> DiffusionStrideCases()
    {
        foreach (var property in typeof(ErrorDiffusionDithering).GetProperties(System.Reflection.BindingFlags.Static | System.Reflection.BindingFlags.Public))
            foreach (bool negative in new[] { false, true }) yield return new object[] { property.Name, negative };
    }

    [Theory]
    [MemberData(nameof(DiffusionStrideCases))]
    public void EveryDiffusionKernelIsIndependentOfPhysicalRowLayout(string name, bool negative)
    {
        const int width = 7, height = 5, pitch = 48;
        var property = typeof(ErrorDiffusionDithering).GetProperty(name)!;
        ErrorDiffusionDithering Filter() => (ErrorDiffusionDithering)property.GetValue(null)!;
        using var original = CreateImage(width, height, (x, y) => Color.FromArgb(11 + x + y, (x * 41 + y * 13) % 256, (x * 29 + y * 17) % 256, (x * 19 + y * 53) % 256));
        using var expected = (Bitmap)original.Clone(); Filter().Apply(expected);
        var bytes = Enumerable.Repeat((byte)173, 3 * height * pitch).ToArray();
        var active = new bool[bytes.Length];
        int step = negative ? -pitch : pitch, origin = height * pitch + (negative ? (height - 1) * pitch : 0);
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++)
        {
            var color = original.GetPixel(x, y);
            byte[] pixel = { color.B, color.G, color.R, color.A };
            for (int c = 0; c < 4; c++) { int i = origin + y * step + x * 4 + c; bytes[i] = pixel[c]; active[i] = true; }
        }
        var handle = GCHandle.Alloc(bytes, GCHandleType.Pinned);
        try
        {
            Filter().Apply(new BitmapData { Width = width, Height = height, Stride = step,
                PixelFormat = PixelFormat.Format32bppArgb, Scan0 = IntPtr.Add(handle.AddrOfPinnedObject(), origin) });
        }
        finally { handle.Free(); }
        for (int i = 0; i < bytes.Length; i++) if (!active[i]) Assert.Equal((byte)173, bytes[i]);
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++)
        {
            int i = origin + y * step + x * 4;
            ImagingAuditTests.Pixel(expected.GetPixel(x, y), Color.FromArgb(bytes[i + 3], bytes[i + 2], bytes[i + 1], bytes[i]));
        }
    }

    /// <summary>Applies a scalar Floyd-Steinberg raster to B,G,R,A samples, preserving alpha</summary>
    /// <param name="pixels">Integer channels, updated in place with truncation after each diffusion step.</param>
    /// <param name="width">Logical width.</param>
    /// <param name="height">Logical height.</param>
    private static void ReferenceFloydSteinberg(int[,,] pixels, int width, int height)
    {
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++) for (int c = 0; c < 3; c++)
        {
            int old = pixels[y, x, c], value = old < 129 ? 0 : 255;
            pixels[y, x, c] = value;
            foreach (var (dx, dy, weight) in new[] { (1, 0, 7), (-1, 1, 3), (0, 1, 5), (1, 1, 1) })
            {
                int xx = x + dx, yy = y + dy;
                if (xx >= 0 && xx < width && yy < height)
                    pixels[yy, xx, c] = QuantizeByte(pixels[yy, xx, c] + (old - value) * weight / 16.0);
            }
        }
    }

    [Theory]
    [InlineData(0)] [InlineData(1)] [InlineData(7)] [InlineData(41)] [InlineData(103)]
    public void HistogramMediansMatchSortedObservationsForOddAndEvenPopulations(int seed)
    {
        var random = new Random(seed);
        for (int count = 0; count <= 129; count++)
        {
            int[] samples = Enumerable.Range(0, count).Select(_ => random.Next(256)).OrderBy(x => x).ToArray();
            var histogram = new int[256];
            foreach (int value in samples) histogram[value]++;
            Assert.Equal(count == 0 ? 0 : samples[(count - 1) / 2], Statistics.Median(histogram));
        }
    }

    [Theory]
    [InlineData(0)] [InlineData(20)] [InlineData(255)]
    public void ARepeatedSingleIntensityIsItsOwnMedian(int value)
    {
        foreach (int count in new[] { 1, 2, 3, 10, int.MaxValue })
        {
            var histogram = new int[256]; histogram[value] = count;
            Assert.Equal(value, Statistics.Median(histogram));
        }
    }

    [Fact]
    public void HistogramRanksRemainExactBeyondInt32PopulationLimits()
    {
        Assert.Equal(2, Statistics.Median(new[] { 0, int.MaxValue, int.MaxValue, int.MaxValue }));
        Assert.Equal(1, Statistics.Median(new[] { 0, int.MaxValue, 0, int.MaxValue }));
        Assert.Equal(0, Statistics.Median(Array.Empty<int>()));
        Assert.Equal(0, Statistics.Median(new int[256]));
        Assert.Throws<ArgumentException>(() => Statistics.Median(new[] { 2, -1, 3 }));
    }

    public static IEnumerable<object[]> TransferCases()
    {
        foreach (string target in new[] { "constant", "horizontal", "vertical", "mixed", "constant-channel" })
            foreach (string source in new[] { "constant", "horizontal", "vertical", "mixed", "constant-channel" })
                foreach (bool inverted in new[] { false, true }) foreach (float factor in new[] { 0f, .5f, 2f })
                    yield return new object[] { target, source, inverted, factor };
    }

    [Theory]
    [MemberData(nameof(TransferCases))]
    public void ColorTransferMatchesIndependentPopulationMoments(string targetKind, string sourceKind, bool inverted, float factor)
    {
        using var target = CreateImage(5, 3, (x, y) => TransferPixel(x, y, targetKind, 0));
        using var source = CreateImage(3, 4, (x, y) => TransferPixel(x, y, sourceKind, 1));
        using var sourceCopy = (Bitmap)source.Clone();
        using var expected = ReferenceTransfer(target, source, inverted, factor);
        new ColorTransfer(factor, inverted).Apply(target, source);
        ImagingAuditTests.Same(expected, target, 1);
        ImagingAuditTests.Same(sourceCopy, source);
    }

    public static IEnumerable<object[]> TransferIdentityCases()
    {
        foreach (Space space in new[] { Space.RGB, Space.HSB, Space.HSL, Space.YCbCr })
            foreach (bool inverted in new[] { false, true }) foreach (bool constant in new[] { false, true })
                foreach (var size in new[] { (1, 1), (1, 9), (9, 1), (7, 5) })
                    yield return new object[] { space, inverted, constant, size.Item1, size.Item2 };
    }

    [Theory]
    [MemberData(nameof(TransferIdentityCases))]
    public void SameImageTransferPreservesSingletonThinConstantAndTexturedImages(Space space, bool inverted, bool constant, int width, int height)
    {
        using var expected = CreateImage(width, height, (x, y) => constant ? Color.FromArgb(77, 99, 123) : TransferPixel(x, y, "mixed", 0));
        using var actual = (Bitmap)expected.Clone();
        new ColorTransfer(0, inverted, space).Apply(actual, expected);
        ImagingAuditTests.Same(expected, actual, space is Space.HSB or Space.HSL ? 5 : 1);
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void TransferDependsOnPixelDistributionRatherThanReferenceShapeOrReplication(bool inverted)
    {
        using var original = CreateImage(5, 3, (x, y) => TransferPixel(x, y, "mixed", 0));
        using var source = CreateImage(3, 5, (x, y) => TransferPixel(x, y, "mixed", 1));
        using var reshaped = CreateImage(5, 3, (x, y) => source.GetPixel((x + 5 * y) % 3, (x + 5 * y) / 3));
        using var repeated = CreateImage(6, 10, (x, y) => source.GetPixel(x % 3, y % 5));
        using var expected = (Bitmap)original.Clone(); new ColorTransfer(0, inverted).Apply(expected, source);
        foreach (var reference in new[] { reshaped, repeated })
        {
            using var actual = (Bitmap)original.Clone(); new ColorTransfer(0, inverted).Apply(actual, reference);
            ImagingAuditTests.Same(expected, actual, 1);
        }
    }

    /// <summary>Builds exact byte samples with varied spatial and per-channel variance</summary>
    /// <param name="x">Column coordinate.</param><param name="y">Row coordinate.</param>
    /// <param name="kind">The spatial pattern, including constant-channel degeneracies.</param>
    /// <param name="reference">Selects a distinct mean and contrast for the reference image.</param>
    /// <returns>An opaque test pixel.</returns>
    private static Color TransferPixel(int x, int y, string kind, int reference)
    {
        int v = kind switch { "constant" => 0, "horizontal" => x * 9, "vertical" => y * 11, _ => (x * 17 + y * 23) % 61 };
        int r = kind == "constant-channel" ? 113 : 61 + v;
        return Color.FromArgb(r + reference * 27, 79 + v / 2 + reference * 11, 93 + v + reference * 17);
    }

    /// <summary>Computes a byte-space transfer reference using pairwise population variance</summary>
    /// <param name="target">Destination samples before transfer.</param><param name="source">Reference samples.</param>
    /// <param name="inverted">Whether to use the existing reciprocal gain.</param><param name="factor">Contrast factor.</param>
    /// <returns>An independently computed opaque bitmap, saturated and truncated to bytes.</returns>
    private static Bitmap ReferenceTransfer(Bitmap target, Bitmap source, bool inverted, float factor)
    {
        double[][] Samples(Bitmap image) => new Func<Color, int>[] { c => c.R, c => c.G, c => c.B }
            .Select(component => Enumerable.Range(0, image.Width * image.Height)
                .Select(i => (double)component(image.GetPixel(i % image.Width, i / image.Width))).ToArray()).ToArray();
        double Deviation(double[] values)
        {
            // Var(X) = sum(i<j, (x_i-x_j)^2) / N^2, independent of Welford's recurrence.
            double sum = 0;
            for (int i = 0; i < values.Length; i++) for (int j = i + 1; j < values.Length; j++)
                sum += (values[i] - values[j]) * (values[i] - values[j]);
            return Math.Sqrt(sum) / values.Length;
        }
        var a = Samples(target); var b = Samples(source);
        var output = new int[3][];
        for (int c = 0; c < 3; c++)
        {
            double da = Deviation(a[c]), db = Deviation(b[c]), mean = a[c].Average(), referenceMean = b[c].Average();
            double gain = da == 0 || db == 0 ? 0 : inverted ? db / da / (1 + factor) : da / db * (1 + factor);
            output[c] = a[c].Select(v => QuantizeByte(referenceMean + gain * (v - mean))).ToArray();
        }
        return CreateImage(target.Width, target.Height, (x, y) =>
        {
            int i = y * target.Width + x;
            return Color.FromArgb(output[0][i], output[1][i], output[2][i]);
        });
    }

    /// <summary>Clamps a finite scalar to the byte range and truncates its fractional part</summary>
    /// <param name="value">The unquantized sample.</param><returns>An integer from zero through 255.</returns>
    private static int QuantizeByte(double value) => (int)Math.Clamp(value, 0, 255);

    /// <summary>Creates a 32-bit test image from independently specified pixel values</summary>
    /// <param name="width">Positive width.</param><param name="height">Positive height.</param>
    /// <param name="pixel">Pixel generator addressed by column and row.</param><returns>The new bitmap.</returns>
    private static Bitmap CreateImage(int width, int height, Func<int, int, Color> pixel)
    {
        var bitmap = new Bitmap(width, height, PixelFormat.Format32bppArgb);
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++) bitmap.SetPixel(x, y, pixel(x, y));
        return bitmap;
    }
}
