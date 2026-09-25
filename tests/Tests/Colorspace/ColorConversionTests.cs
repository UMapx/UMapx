using System.Reflection;
using UMapx.Colorspace;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Colorspace")]
public class ColorConversionTests
{
    private static readonly string[] Models =
    {
        "sRGB",
        "CMYK",
        "HSB",
        "HSL",
        "AHSL",
        "RYB",
        "XYZ",
        "LAB",
        "YUV",
        "YPbPr",
        "YIQ",
        "YDbDr",
        "YCgCo",
        "YCbCr"
    };
    public static IEnumerable<object[]> ColorCases()
    {
        foreach (string name in Models)
            foreach (int r in new[]
            {
                0,
                64,
                128,
                255
            }

            )
                foreach (int g in new[]
                {
                    0,
                    64,
                    128,
                    255
                }

                )
                    foreach (int b in new[]
                    {
                        0,
                        64,
                        128,
                        255
                    }

                    )
                        yield return new object[]
                        {
                            name,
                            r,
                            g,
                            b
                        };
    }

    private static void CloseColor(RGB expected, RGB actual, int tolerance)
    {
        Assert.True(Math.Abs(expected.Red - actual.Red) <= tolerance && Math.Abs(expected.Green - actual.Green) <= tolerance && Math.Abs(expected.Blue - actual.Blue) <= tolerance, $"Expected RGB({expected.Red}, {expected.Green}, {expected.Blue}); actual RGB({actual.Red}, {actual.Green}, {actual.Blue}); channel tolerance {tolerance}.");
    }

    [Theory, MemberData(nameof(ColorCases))]
    public void ColorConversionsPreserveColorsWithinTheirQuantizationBudget(string name, int r, int g, int b)
    {
        var type = typeof(RGB).Assembly.GetType("UMapx.Colorspace." + name)!;
        var rgb = new RGB(r, g, b);
        var converted = (IColorSpace)type.GetMethod("FromRGB", new[] { typeof(int), typeof(int), typeof(int) })!.Invoke(null, new object[] { r, g, b })!;
        var overload = type.GetMethod("FromRGB", new[] { typeof(RGB) })!.Invoke(null, new object[] { rgb })!;
        Assert.Equal(converted, overload);
        // HSB/HSL intentionally quantize hue to whole degrees: at most 255/60 + one output byte.
        CloseColor(rgb, converted.ToRGB, name is "HSB" or "HSL" ? 6 : 2);
    }

    [Theory]
    [InlineData(0, 0, 0)]
    [InlineData(255, 255, 255)]
    [InlineData(255, 0, 0)]
    [InlineData(0, 255, 0)]
    [InlineData(0, 0, 255)]
    [InlineData(12, 64, 190)]
    [InlineData(128, 0, 255)]
    public void StandardColorComponentsAgreeWithIndependentFormulas(int red, int green, int blue)
    {
        double r = red / 255.0, g = green / 255.0, b = blue / 255.0, max = Math.Max(r, Math.Max(g, b)), min = Math.Min(r, Math.Min(g, b)), delta = max - min;
        double hue = delta == 0 ? 0 : max == r ? 60 * ((g - b) / delta + 6) % 360 : max == g ? 60 * ((b - r) / delta + 2) : 60 * ((r - g) / delta + 4);
        var hsv = HSB.FromRGB(red, green, blue);
        Close(hue, hsv.Hue, 1.0001, 0);
        Close(max, hsv.Brightness);
        Close(max == 0 ? 0 : delta / max, hsv.Saturation);
        var hsl = HSL.FromRGB(red, green, blue);
        double light = (max + min) / 2;
        Close(hue, hsl.Hue, 1.0001, 0);
        Close(light, hsl.Lightness);
        Close(delta == 0 ? 0 : delta / (1 - Math.Abs(2 * light - 1)), hsl.Saturation);
        var cmyk = CMYK.FromRGB(red, green, blue);
        Close(1 - max, cmyk.Keycolor);
        Close(max == 0 ? 0 : 1 - r / max, cmyk.Cyan);
        Close(max == 0 ? 0 : 1 - g / max, cmyk.Magenta);
        Close(max == 0 ? 0 : 1 - b / max, cmyk.Yellow);
        var srgb = sRGB.FromRGB(red, green, blue);
        Close(r, srgb.Red);
        Close(g, srgb.Green);
        Close(b, srgb.Blue);
        double y = .299 * r + .587 * g + .114 * b;
        var cbcr = YCbCr.FromRGB(red, green, blue);
        Close(y, cbcr.Y);
        Close(.5 + (b - y) / (2 * (1 - .114)), cbcr.Cb, 2e-6);
        Close(.5 + (r - y) / (2 * (1 - .299)), cbcr.Cr, 2e-6);
        var cgco = YCgCo.FromRGB(red, green, blue);
        Close((r + 2 * g + b) / 4, cgco.Y);
        Close((-r + 2 * g - b) / 4, cgco.Cg);
        Close((r - b) / 2, cgco.Co);
    }

    [Theory]
    [InlineData(0, 0, 0)]
    [InlineData(255, 255, 255)]
    [InlineData(255, 0, 0)]
    [InlineData(0, 255, 0)]
    [InlineData(0, 0, 255)]
    [InlineData(12, 64, 190)]
    [InlineData(128, 128, 128)]
    public void CieXyzAndLabMatchTheD65ReferenceWhite(int red, int green, int blue)
    {
        double Linear(int v) => v / 255.0 <= .04045 ? v / 255.0 / 12.92 : Math.Pow((v / 255.0 + .055) / 1.055, 2.4);
        double r = Linear(red), g = Linear(green), b = Linear(blue);
        double x = .4124 * r + .3576 * g + .1805 * b, y = .2126 * r + .7152 * g + .0722 * b, z = .0193 * r + .1192 * g + .9505 * b;
        var xyz = XYZ.FromRGB(red, green, blue);
        Close(x, xyz.X, 1e-5);
        Close(y, xyz.Y, 1e-5);
        Close(z, xyz.Z, 1e-5);
        double F(double t) => t > Math.Pow(6.0 / 29, 3) ? Math.Cbrt(t) : t / (3 * Math.Pow(6.0 / 29, 2)) + 4.0 / 29;
        double l = 116 * F(y) - 16, a = 500 * (F(x / .9505) - F(y)), bb = 200 * (F(y) - F(z / 1.089));
        var lab = LAB.FromRGB(red, green, blue);
        Close(l, lab.L, .003);
        Close(a, lab.A, .003);
        Close(bb, lab.B, .003);
        var restored = LAB.ToXYZ((float)l, (float)a, (float)bb);
        Close(x, restored.X, .0001);
        Close(y, restored.Y, .0001);
        Close(z, restored.Z, .0001);
    }

    public static IEnumerable<object[]> ModelCases() => Models.Select(n => new object[] { n });
    [Theory, MemberData(nameof(ModelCases))]
    public void ColorValuesCloneAndPreserveTheirComponentsInUnknownContainers(string name)
    {
        var type = typeof(RGB).Assembly.GetType("UMapx.Colorspace." + name)!;
        object value = type.GetMethod("FromRGB", new[] { typeof(RGB) })!.Invoke(null, new object[] { new RGB(80, 120, 170) })!;
        var clone = type.GetMethod("Clone", Type.EmptyTypes)!.Invoke(value, null)!;
        Assert.Equal(value, clone);
        Assert.Equal(value.GetHashCode(), clone.GetHashCode());
        Assert.Equal(value, ((ICloneable)value).Clone());
        Assert.False(value.Equals(null));
        Assert.False(value.Equals("color"));
        Assert.False(string.IsNullOrWhiteSpace(value.ToString()));
        Assert.True((bool)type.GetMethod("op_Equality")!.Invoke(null, new[] { value, clone })!);
        Assert.False((bool)type.GetMethod("op_Inequality")!.Invoke(null, new[] { value, clone })!);
        if (name == "CMYK")
            return; // A three-component container cannot preserve a fourth independent channel.
        var to = typeof(Unknown).GetMethods(BindingFlags.Public | BindingFlags.Static).Single(m => m.Name == "op_Implicit" && m.ReturnType == typeof(Unknown) && m.GetParameters()[0].ParameterType == type);
        var from = typeof(Unknown).GetMethods(BindingFlags.Public | BindingFlags.Static).Single(m => m.Name == "op_Implicit" && m.ReturnType == type);
        object container = to.Invoke(null, new[] { value })!;
        Assert.Equal(value, from.Invoke(null, new[] { container }));
        var u = (Unknown)container;
        Assert.Equal(u, u.Clone());
        Assert.True(u == u.Clone());
        Assert.False(u != u.Clone());
        Assert.Equal(u.GetHashCode(), u.Clone().GetHashCode());
    }

    [Theory]
    [InlineData(0, 0, 0)]
    [InlineData(255, 255, 255)]
    [InlineData(12, 34, 56)]
    [InlineData(128, 1, 254)]
    public void RgbTextConversionsAndNeutralAdjustmentsPreserveBytes(int r, int g, int b)
    {
        var rgb = new RGB(r, g, b);
        var hex = RGB.ToHEX(rgb);
        Assert.Equal(rgb, RGB.FromHEX(hex));
        Assert.Equal(hex, RGB.ToHEX(r, g, b));
        Assert.Equal(rgb, RGB.FromHEX($"#{r:X2}{g:X2}{b:X2}"));
        CloseColor(rgb, RGB.Saturation(rgb, 0), 1);
        CloseColor(rgb, RGB.Vibrance(rgb, 0), 1);
        Close((r + g + b) / 3, RGB.Average(rgb));
        Close((r + g + b) / 3.0 / 255, RGB.Average(sRGB.FromRGB(rgb)));
    }

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
            "XYZ" => XYZ.FromRGB(value),
            "LAB" => LAB.FromRGB(value),
            "RYB" => RYB.FromRGB(value),
            _ => AHSL.FromRGB(value)
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
        Assert.Equal(x, value.X);
        Assert.Equal(y, value.Y);
        Assert.Equal(z, value.Z);
        Assert.Equal(value, value.Clone());
        var assigned = new XYZ
        {
            X = x,
            Y = y,
            Z = z
        };
        Assert.Equal(value, assigned);
        Assert.Equal(new XYZ(), new XYZ(-1, -2, -3));
    }

    [Theory]
    [InlineData(0f)]
    [InlineData(1f)]
    [InlineData(7.999f)]
    [InlineData(8f)]
    [InlineData(8.001f)]
    [InlineData(50f)]
    [InlineData(100f)]
    public void LabInverseIsContinuousAcrossThePiecewiseCieThreshold(float lightness)
    {
        var xyz = LAB.ToXYZ(lightness, 0, 0);
        double y = lightness > 8 ? Math.Pow((lightness + 16.0) / 116, 3) : lightness * 27.0 / 24389;
        Close(.9505 * y, xyz.X, 3e-7, 0);
        Close(y, xyz.Y, 3e-7, 0);
        Close(1.089 * y, xyz.Z, 3e-7, 0);
        var lab = XYZ.ToLAB(xyz);
        Close(lightness, lab.L, 2e-5, 0);
        Close(0, lab.A, 6e-5, 0);
        Close(0, lab.B, 6e-5, 0);
    }

    private static void Check(RGB expected, RGB actual, int tolerance)
    {
        Assert.True(Math.Abs(expected.Red - actual.Red) <= tolerance && Math.Abs(expected.Green - actual.Green) <= tolerance && Math.Abs(expected.Blue - actual.Blue) <= tolerance, $"RGB({expected.Red},{expected.Green},{expected.Blue}) -> RGB({actual.Red},{actual.Green},{actual.Blue})");
    }

    [Fact]
    public void XyzCanRepresentItsD65WhitePoint() => Close(1.089, XYZ.White.Z);
    [Fact]
    public void RybConversionPreservesWhite()
    {
        var color = RYB.FromRGB(255, 255, 255).ToRGB;
        Assert.Equal((byte)255, color.Red);
        Assert.Equal((byte)255, color.Green);
        Assert.Equal((byte)255, color.Blue);
    }
}
