using UMapx.Colorspace;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "ColorSpace")]
public class ColorSpaceRepairTests
{
    [Fact]
    public void EveryNeutralBytePreservesWhiteAndHasZeroLabChroma()
    {
        for (int v = 0; v <= 255; v++)
        {
            var input = new RGB(v, v, v);
            Check(input, RYB.FromRGB(input).ToRGB, 1);
            Check(input, AHSL.FromRGB(input).ToRGB, 1);
            var lab = LAB.FromRGB(input);
            Close(0, lab.A, 6e-5, 0);
            Close(0, lab.B, 6e-5, 0);
            Check(input, lab.ToRGB, 1);
        }
    }

    [Theory]
    [InlineData("XYZ")]
    [InlineData("LAB")]
    [InlineData("RYB")]
    [InlineData("AHSL")]
    public void DeterministicColorsAndRedEqualToGrayKeepTheirChroma(string model)
    {
        IColorSpace Convert(RGB value) => model switch
        {
            "XYZ" => XYZ.FromRGB(value), "LAB" => LAB.FromRGB(value),
            "RYB" => RYB.FromRGB(value), _ => AHSL.FromRGB(value)
        };
        var random = new Random(7010);
        for (int i = 0; i < 1024; i++)
        {
            var input = new RGB(random.Next(256), random.Next(256), random.Next(256));
            Check(input, Convert(input).ToRGB, 2);
        }
        for (int red = 1; red < 255; red += 7)
        for (int delta = 1; delta <= Math.Min(red, 255 - red); delta += 5)
        {
            var input = new RGB(red, red - delta, red + delta);
            Check(input, Convert(input).ToRGB, 2);
            Assert.True(AHSL.FromRGB(input).Saturation > 0);
        }
    }

    [Theory]
    [InlineData(.9505f, 1f, 1.089f)]
    [InlineData(1.4f, 2f, 2.8f)]
    [InlineData(0f, 0f, 0f)]
    public void XyzCoordinatesAreNonnegativeAndRelativeRatherThanClippedToOne(float x, float y, float z)
    {
        var value = new XYZ(x, y, z);
        Assert.Equal(x, value.X); Assert.Equal(y, value.Y); Assert.Equal(z, value.Z);
        Assert.Equal(value, value.Clone());
        var assigned = new XYZ { X = x, Y = y, Z = z };
        Assert.Equal(value, assigned);
        Assert.Equal(new XYZ(), new XYZ(-1, -2, -3));
    }

    [Theory]
    [InlineData(0f)] [InlineData(1f)] [InlineData(7.999f)]
    [InlineData(8f)] [InlineData(8.001f)] [InlineData(50f)] [InlineData(100f)]
    public void LabInverseIsContinuousAcrossThePiecewiseCieThreshold(float lightness)
    {
        var xyz = LAB.ToXYZ(lightness, 0, 0);
        double y = lightness > 8 ? Math.Pow((lightness + 16.0) / 116, 3) : lightness * 27.0 / 24389;
        Close(.9505 * y, xyz.X, 3e-7, 0);
        Close(y, xyz.Y, 3e-7, 0);
        Close(1.089 * y, xyz.Z, 3e-7, 0);
        var lab = XYZ.ToLAB(xyz);
        Close(lightness, lab.L, 2e-5, 0);
        Close(0, lab.A, 6e-5, 0); Close(0, lab.B, 6e-5, 0);
    }

    private static void Check(RGB expected, RGB actual, int tolerance)
    {
        Assert.True(Math.Abs(expected.Red - actual.Red) <= tolerance &&
            Math.Abs(expected.Green - actual.Green) <= tolerance && Math.Abs(expected.Blue - actual.Blue) <= tolerance,
            $"RGB({expected.Red},{expected.Green},{expected.Blue}) -> RGB({actual.Red},{actual.Green},{actual.Blue})");
    }
}
