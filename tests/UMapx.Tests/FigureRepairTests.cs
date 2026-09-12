using System.Drawing;
using System.Runtime.Versioning;
using UMapx.Core;
using UMapx.Visualization;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Geometry")]
[SupportedOSPlatform("windows")]
public class FigureRepairTests
{
    public static IEnumerable<object[]> ConstantCases()
    {
        foreach (var series in Enum.GetValues<SeriesType>())
            foreach (string axis in new[] { "x", "y", "both" })
                foreach (float value in new[] { 0f, 2f, -3f, float.Epsilon, -float.Epsilon, 1e20f, -1e20f, float.MaxValue, -float.MaxValue })
                    yield return new object[] { series, axis, value };
    }

    [Theory]
    [MemberData(nameof(ConstantCases))]
    public void ConstantAxesRenderFiniteRangesAndVisibleMarkers(SeriesType series, string axis, float value)
    {
        using var style = FigureStyle.Standard;
        var figure = new Figure(style); figure.Grid.Show = true; figure.Legend.Show = false;
        float[] varying = { -1, 0, 1 };
        float[] constant = { value, value, value };
        var x = axis == "y" ? varying : constant;
        var y = axis == "x" ? varying : constant;
        figure.Plot(new PlotSeries(x, y, 2, Color.Red, series, ShapeType.Circle, "constant"));
        using var bitmap = new Bitmap(480, 320); figure.To(bitmap);
        AssertRange(figure.RangeX, x); AssertRange(figure.RangeY, y);
        if (axis == "x") Assert.Equal((-1f, 1f), (figure.RangeY.Min, figure.RangeY.Max));
        if (axis == "y") Assert.Equal((-1f, 1f), (figure.RangeX.Min, figure.RangeX.Max));
        Assert.True(CountRed(bitmap) > 5, "The constant series has no visible markers.");
        var rangeX = figure.RangeX; var rangeY = figure.RangeY;
        figure.To(bitmap);
        Assert.Equal((rangeX.Min, rangeX.Max), (figure.RangeX.Min, figure.RangeX.Max));
        Assert.Equal((rangeY.Min, rangeY.Max), (figure.RangeY.Min, figure.RangeY.Max));
    }

    [Theory]
    [InlineData(-3f)] [InlineData(0f)] [InlineData(2f)]
    public void SingletonSeriesCanBeRenderedWithoutManualAxisLimits(float value)
    {
        using var style = FigureStyle.Standard; var figure = new Figure(style); figure.Legend.Show = false;
        figure.Plot(new PlotSeries(new[] { value }, new[] { value }, 2, Color.Red, SeriesType.Scatter, ShapeType.Ball, "single"));
        using var bitmap = new Bitmap(320, 240); figure.To(bitmap);
        AssertRange(figure.RangeX, new[] { value }); AssertRange(figure.RangeY, new[] { value });
        Assert.True(CountRed(bitmap) > 5);
    }

    [Fact]
    public void MultipleConstantSeriesUseTheirCombinedDataBounds()
    {
        using var style = FigureStyle.Standard; var figure = new Figure(style);
        figure.Plot(new PlotSeries(new[] { 1f, 2f }, new[] { -3f, -3f }, 2, Color.Red, SeriesType.Plot, ShapeType.None, "lower"));
        figure.Plot(new PlotSeries(new[] { 4f, 6f }, new[] { 7f, 7f }, 2, Color.Blue, SeriesType.Plot, ShapeType.None, "upper"));
        using var bitmap = new Bitmap(480, 320); figure.To(bitmap);
        Assert.Equal((1f, 6f), (figure.RangeX.Min, figure.RangeX.Max));
        Assert.Equal((-3f, 7f), (figure.RangeY.Min, figure.RangeY.Max));
    }

    [Fact]
    public void ConstantSeriesDoNotOverrideManualRanges()
    {
        using var style = FigureStyle.Standard;
        var figure = new Figure(style) { AutoRange = false, RangeX = new RangeFloat(-10, 10), RangeY = new RangeFloat(-5, 5) };
        figure.Plot(new PlotSeries(new[] { 2f, 2f }, new[] { 2f, 2f }, 2, Color.Red, SeriesType.Scatter, ShapeType.Ball, "constant"));
        using var bitmap = new Bitmap(320, 240); figure.To(bitmap);
        Assert.Equal((-10f, 10f), (figure.RangeX.Min, figure.RangeX.Max));
        Assert.Equal((-5f, 5f), (figure.RangeY.Min, figure.RangeY.Max));
        Assert.Throws<ArgumentOutOfRangeException>(() => figure.RangeX = new RangeFloat(2, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => figure.RangeY = new RangeFloat(2, 2));
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void EmptyOrNonfiniteDataRetainUsableFallbackRanges(bool nonfinite)
    {
        using var style = FigureStyle.Standard; var figure = new Figure(style);
        var values = nonfinite ? new[] { float.NaN, float.PositiveInfinity, float.NegativeInfinity } : Array.Empty<float>();
        figure.Plot(new PlotSeries(values, values, 2, Color.Red, SeriesType.Plot, ShapeType.None, "empty"));
        using var bitmap = new Bitmap(320, 240); figure.To(bitmap);
        Assert.Equal((-5f, 5f), (figure.RangeX.Min, figure.RangeX.Max));
        Assert.Equal((-5f, 5f), (figure.RangeY.Min, figure.RangeY.Max));
    }

    /// <summary>Checks finite, nonzero axis extent and inclusion of every supplied finite sample</summary>
    /// <param name="range">Computed axis range.</param><param name="samples">Finite data values.</param>
    private static void AssertRange(RangeFloat range, float[] samples)
    {
        Assert.True(float.IsFinite(range.Min) && float.IsFinite(range.Max));
        Assert.True(range.Min < range.Max);
        foreach (float value in samples) Assert.InRange(value, range.Min, range.Max);
    }

    /// <summary>Counts series-colored pixels independently of the figure's coordinate mapping</summary>
    /// <param name="bitmap">Rendered figure with its legend disabled.</param><returns>The number of red pixels.</returns>
    private static int CountRed(Bitmap bitmap)
    {
        int count = 0;
        for (int y = 0; y < bitmap.Height; y++) for (int x = 0; x < bitmap.Width; x++)
        {
            var color = bitmap.GetPixel(x, y);
            if (color.R > 180 && color.G < 90 && color.B < 90) count++;
        }
        return count;
    }
}
