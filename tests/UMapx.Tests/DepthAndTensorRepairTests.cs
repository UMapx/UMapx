using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.Versioning;
using UMapx.Imaging;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Imaging")]
[SupportedOSPlatform("windows")]
public class DepthAndTensorRepairTests
{
    public static IEnumerable<object[]> TensorCases()
    {
        foreach (bool rgb in new[] { false, true }) foreach (bool floating in new[] { false, true })
            foreach (var size in new[] { (1, 1), (1, 9), (9, 1), (7, 5), (256, 1) })
                yield return new object[] { rgb, floating, size.Item1, size.Item2 };
    }

    [Theory]
    [MemberData(nameof(TensorCases))]
    public void ThreePlaneTensorsProduceOpaquePixelsWithTheRequestedChannelOrder(bool rgb, bool floating, int width, int height)
    {
        var bytes = Enumerable.Range(0, 3).Select(c => Enumerable.Range(0, width * height)
            .Select(i => (byte)((i * 73 + c * 59) % 256)).ToArray()).ToArray();
        var values = bytes.Select(plane => plane.Select(v => Math.Min(255f, v + .75f)).ToArray()).ToArray();
        using var image = floating ? values.FromFloatTensor(width, height, rgb) : bytes.FromByteTensor(width, height, rgb);
        Assert.Equal(new Size(width, height), image.Size);
        Assert.Equal(PixelFormat.Format32bppArgb, image.PixelFormat);
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++)
        {
            int i = y * width + x;
            var expected = Color.FromArgb(255, bytes[rgb ? 0 : 2][i], bytes[1][i], bytes[rgb ? 2 : 0][i]);
            ImagingAuditTests.Pixel(expected, image.GetPixel(x, y));
        }
    }

    public static IEnumerable<object[]> ConstantDepthCases()
    {
        foreach (int value in new[] { 0, 1, 1024, 65535 })
            foreach (var size in new[] { (1, 1), (65535, 1), (256, 256), (300, 300), (257, 511) })
                yield return new object[] { value, size.Item1, size.Item2 };
    }

    [Theory]
    [MemberData(nameof(ConstantDepthCases))]
    public void ConstantDepthCdfReachesTheMaximumAcrossCountOverflowBoundaries(int value, int width, int height)
    {
        var input = CreateDepth(width, height, (x, y) => (ushort)value);
        var result = input.Equalize();
        Assert.Equal(height, result.GetLength(0)); Assert.Equal(width, result.GetLength(1));
        foreach (ushort output in result) Assert.Equal(ushort.MaxValue, output);
        foreach (ushort original in input) Assert.Equal((ushort)value, original);
    }

    public static IEnumerable<object[]> DepthDistributionCases()
    {
        foreach (int kind in new[] { 0, 1, 2 })
            foreach (var size in new[] { (7, 5), (257, 263), (65537, 1), (1, 65537) })
                yield return new object[] { kind, size.Item1, size.Item2 };
    }

    [Theory]
    [MemberData(nameof(DepthDistributionCases))]
    public void DepthEqualizationMatchesExactRanksOfSortedObservations(int kind, int width, int height)
    {
        var input = CreateDepth(width, height, (x, y) =>
        {
            int i = y * width + x;
            return kind switch { 0 => (ushort)((i * 251L) % 65536), 1 => (ushort)(i % 3 * 32767), _ => i % 31 == 0 ? ushort.MaxValue : (ushort)17 };
        });
        var sorted = input.Cast<ushort>().OrderBy(v => v).ToArray();
        var original = (ushort[,])input.Clone(); var actual = input.Equalize();
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++)
        {
            int rank = UpperBound(sorted, input[y, x]);
            Assert.Equal((ushort)(65535L * rank / sorted.Length), actual[y, x]);
            Assert.Equal(original[y, x], input[y, x]);
        }
    }

    [Theory]
    [InlineData(1)] [InlineData(2)]
    public void EveryDepthCodeMatchesItsExactUniformCdf(int repetitions)
    {
        var input = CreateDepth(65536, repetitions, (x, y) => (ushort)x);
        var output = input.Equalize();
        for (int y = 0; y < repetitions; y++) for (int x = 0; x < 65536; x++)
            Assert.Equal((ushort)((x + 1L) * 65535 / 65536), output[y, x]);
    }

    [Theory]
    [InlineData(0, 0)] [InlineData(0, 7)] [InlineData(7, 0)]
    public void EmptyDepthEqualizationPreservesShape(int width, int height)
    {
        var output = new ushort[height, width].Equalize();
        Assert.Equal(height, output.GetLength(0)); Assert.Equal(width, output.GetLength(1));
    }

    public static IEnumerable<object[]> MergeCases()
    {
        foreach (var size in new[] { (1, 1), (4, 3), (3, 4) })
            foreach (var origin in new[] { (0, 0), (2, 1), (-1, 2), (3, -1), (-2, -1), (6, 4), (7, 5), (int.MinValue, 0), (0, int.MinValue), (int.MaxValue, 1), (1, int.MaxValue) })
                yield return new object[] { size.Item1, size.Item2, origin.Item1, origin.Item2 };
    }

    [Theory]
    [MemberData(nameof(MergeCases))]
    public void EqualSizeDepthPlacementsCopyExactSamplesAndClipAllFourEdges(int width, int height, int x, int y)
    {
        var target = CreateDepth(7, 5, (xx, yy) => (ushort)(100 + 11 * xx + 31 * yy));
        var original = (ushort[,])target.Clone();
        var source = CreateDepth(width, height, (xx, yy) => (ushort)(1000 + 101 * xx + 317 * yy));
        var sourceCopy = (ushort[,])source.Clone();
        target.Merge(source, new Rectangle(x, y, width, height));
        AssertPlacement(original, source, target, x, y);
        Assert.Equal(sourceCopy.Cast<ushort>(), source.Cast<ushort>());
    }

    [Theory]
    [InlineData(0, 0)] [InlineData(2, 1)] [InlineData(-2, -1)] [InlineData(6, 4)]
    public void ResizedDepthPlacementClipsAfterTheFullResize(int x, int y)
    {
        var target = CreateDepth(7, 5, (xx, yy) => (ushort)(xx + 10 * yy));
        var original = (ushort[,])target.Clone();
        var source = CreateDepth(3, 2, (xx, yy) => (ushort)(1000 + 100 * xx + 200 * yy));
        // Isolate placement from the existing Resize interpolation convention.
        var resized = source.Resize(new Size(6, 4));
        target.Merge(source, new Rectangle(x, y, 6, 4));
        AssertPlacement(original, resized, target, x, y);
    }

    [Theory]
    [InlineData(1, 0)] [InlineData(0, 1)] [InlineData(-1, -1)]
    public void OverlappingSelfMergesReadTheOriginalDepthSamples(int x, int y)
    {
        var target = CreateDepth(7, 5, (xx, yy) => (ushort)(100 + xx + 10 * yy));
        var original = (ushort[,])target.Clone();
        target.Merge(target, new Rectangle(x, y, 7, 5));
        AssertPlacement(original, original, target, x, y);
    }

    [Fact]
    public void DefaultDepthMergeCopiesSamplesWithoutResampling()
    {
        var target = new ushort[4, 5];
        ushort[,] source = { { 0, 65535, 1 }, { 7, 1000, 17 } };
        var original = (ushort[,])target.Clone();
        target.Merge(source);
        AssertPlacement(original, source, target, 0, 0);
    }

    [Fact]
    public void DepthMergeHandlesEmptyPlacementsAndRejectsNegativeSizes()
    {
        var target = CreateDepth(5, 4, (x, y) => (ushort)(x + 10 * y));
        var original = (ushort[,])target.Clone();
        var source = new ushort[0, 0];
        foreach (var rectangle in new[] { new Rectangle(0, 0, 0, 4), new Rectangle(0, 0, 5, 0), new Rectangle(int.MaxValue, int.MaxValue, 3, 3) })
            target.Merge(source, rectangle);
        new ushort[0, 3].Merge(new ushort[2, 2]);
        Assert.Equal(original.Cast<ushort>(), target.Cast<ushort>());
        Assert.Throws<ArgumentOutOfRangeException>(() => target.Merge(source, new Rectangle(0, 0, -1, 2)));
        Assert.Throws<ArgumentOutOfRangeException>(() => target.Merge(source, new Rectangle(0, 0, 2, -1)));
        Assert.Throws<ArgumentException>(() => target.Merge(source, new Rectangle(0, 0, 2, 2)));
    }

    /// <summary>Checks placement against independent destination-to-source coordinate mapping</summary>
    /// <param name="before">Original destination samples.</param><param name="source">Full placement samples.</param>
    /// <param name="after">Actual destination.</param><param name="x">Placement column.</param><param name="y">Placement row.</param>
    private static void AssertPlacement(ushort[,] before, ushort[,] source, ushort[,] after, int x, int y)
    {
        for (int row = 0; row < before.GetLength(0); row++) for (int col = 0; col < before.GetLength(1); col++)
        {
            long sx = (long)col - x, sy = (long)row - y;
            ushort expected = sx >= 0 && sx < source.GetLength(1) && sy >= 0 && sy < source.GetLength(0)
                ? source[(int)sy, (int)sx] : before[row, col];
            Assert.Equal(expected, after[row, col]);
        }
    }

    /// <summary>Finds the inclusive population rank without using a histogram or the implementation CDF</summary>
    /// <param name="sorted">Ascending samples.</param><param name="value">Observed depth value.</param>
    /// <returns>The number of samples less than or equal to the value.</returns>
    private static int UpperBound(ushort[] sorted, ushort value)
    {
        int first = 0, last = sorted.Length;
        while (first < last)
        {
            int middle = first + (last - first) / 2;
            if (sorted[middle] <= value) first = middle + 1; else last = middle;
        }
        return first;
    }

    /// <summary>Creates a depth map from a coordinate-based sample definition</summary>
    /// <param name="width">Map width.</param><param name="height">Map height.</param>
    /// <param name="sample">Sample generator.</param><returns>The generated map.</returns>
    private static ushort[,] CreateDepth(int width, int height, Func<int, int, ushort> sample)
    {
        var result = new ushort[height, width];
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++) result[y, x] = sample(x, y);
        return result;
    }
}
