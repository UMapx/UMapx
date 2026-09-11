using System.Numerics;
using System.Text.Json;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class ArithmeticRepairTests
{
    public static IEnumerable<object[]> References()
    {
        using var stream = typeof(ArithmeticRepairTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.arithmetic-repair.json")!;
        using var data = JsonDocument.Parse(stream);
        foreach (var row in data.RootElement.GetProperty("cases").EnumerateArray())
            yield return new object[] { row[0].GetString()!, row[1].GetSingle(), row[2].GetSingle(), row[3].GetSingle(), row[4].GetSingle(), row[5].GetDouble(), row[6].GetDouble() };
    }

    [Theory, MemberData(nameof(References))]
    public void ComplexBranchesAndPowersMatchIndependentHighPrecisionReferences(string operation, float re, float im, float exponentRe, float exponentIm, double expectedRe, double expectedIm)
    {
        var z = new Complex32(re, im);
        Complex32 actual = operation switch {
            "Asinh" => Maths.Asinh(z), "Acosh" => Maths.Acosh(z), "Actan" => Maths.Actan(z),
            "Tanh" => Maths.Tanh(z), "Ctanh" => Maths.Ctanh(z), "Sech" => Maths.Sech(z), "Cosch" => Maths.Cosch(z),
            "Pow" => Maths.Pow(re, new Complex32(exponentRe, exponentIm)), _ => throw new ArgumentException(operation) };
        Close(new Complex(expectedRe, expectedIm), actual, 4.0 * float.Epsilon, 3e-6);
    }

    public static IEnumerable<object[]> RealBoundaries()
    {
        foreach (float value in new[] { float.Epsilon, 1e-30f, 1e-10f, 1e-4f, .5f, 1f, 1.0000001f, 20f, 89f, 100f, 1e20f, float.MaxValue })
        {
            yield return new object[] { value };
            yield return new object[] { -value };
        }
    }

    [Theory, MemberData(nameof(RealBoundaries))]
    public void RealHyperbolicFunctionsRemainAccurateAcrossTheFloatRange(float x)
    {
        Rounded(Math.Asinh(x), Maths.Asinh(x));
        Rounded(Math.Tanh(x), Maths.Tanh(x));
        Rounded(1 / Math.Tanh(x), Maths.Ctanh(x));
        Rounded(1 / Math.Sinh(x), Maths.Cosch(x));
        Rounded(1 / Math.Cosh(x), Maths.Sech(x));
        Rounded(Math.Asinh(1.0 / x), Maths.Acosch(x));
        Rounded(Math.Atanh(1.0 / x), Maths.Actanh(x));
        Rounded(Math.Atanh(x), Maths.Atanh(x));
        Rounded(Math.Acosh(x), Maths.Acosh(x));
        Rounded(Math.Acosh(1.0 / x), Maths.Asech(x));
    }

    [Fact]
    public void RealHyperbolicFunctionsPreserveZerosPolesAndDomainErrors()
    {
        float negativeZero = BitConverter.Int32BitsToSingle(int.MinValue);
        foreach (float x in new[] { 0f, negativeZero, float.PositiveInfinity, float.NegativeInfinity, float.NaN })
        {
            Rounded(Math.Asinh(x), Maths.Asinh(x));
            Rounded(Math.Acosh(x), Maths.Acosh(x));
            Rounded(Math.Atanh(x), Maths.Atanh(x));
            Rounded(Math.Tanh(x), Maths.Tanh(x));
            Rounded(1 / Math.Tanh(x), Maths.Ctanh(x));
            Rounded(Math.Asinh(1.0 / x), Maths.Acosch(x));
        }
        Assert.Equal(int.MinValue, BitConverter.SingleToInt32Bits(Maths.Asinh(negativeZero)));
        Assert.Equal(int.MinValue, BitConverter.SingleToInt32Bits(Maths.Atanh(negativeZero)));
        Close(new Complex(Math.PI / 2, 0), Maths.Actan(new Complex32(0, 0)));
    }

    [Theory, MemberData(nameof(RealBoundaries))]
    public void ComplexArithmeticPreservesFiniteResultsAtExtremeScales(float scale)
    {
        var z = new Complex32(scale, -scale);
        Close(Complex.One, z / z, 4.0 * float.Epsilon, 2e-6);
        Close(new Complex(.5, .5), scale / z, 4.0 * float.Epsilon, 2e-6);
        Rounded(Complex.Abs(z), z.Abs);
        Rounded(2.0 * scale * scale, z.Abs2);
        Close(Complex.Log(z), Maths.Log(z), 4.0 * float.Epsilon, 2e-6);
    }

    [Fact]
    public void ComplexProductsAvoidIntermediateOverflowAndPrematureRounding()
    {
        var a = new Complex32(1e30f, 1e30f);
        var b = new Complex32(1e10f, 1e10f);
        var product = a * b;
        Assert.Equal(0, product.Real);
        Assert.Equal(float.PositiveInfinity, product.Imag);
        var zeroDivision = new Complex32(1, 0) / new Complex32(0, 0);
        Assert.True(float.IsNaN(zeroDivision.Real) && float.IsNaN(zeroDivision.Imag));
    }

    public static IEnumerable<object[]> Quadratics()
    {
        foreach (float scale in new[] { 1e-30f, 1e-10f, 1f, 1e10f, 1e30f })
        foreach (var pair in new[] { (0f, 0f), (1f, 1f), (0f, -1f), (0f, 1f), (-3f, 2f), (1e6f, 1f), (-1e6f, 1f) })
            yield return new object[] { scale, pair.Item1 * scale, pair.Item2 * scale };
    }

    [Theory, MemberData(nameof(Quadratics))]
    public void QuadraticRootsSatisfyResidualsAndVietaUnderCoefficientScaling(float a, float b, float c)
    {
        var roots = Maths.Quadratic(a, b, c).Select(z => (Complex)z).ToArray();
        Assert.Equal(2, roots.Length);
        foreach (var root in roots)
        {
            double scale = Math.Abs(a) * root.Magnitude * root.Magnitude + Math.Abs(b) * root.Magnitude + Math.Abs(c);
            Assert.True(Complex.Abs((a * root + b) * root + c) <= 3e-6 * scale + 1e-44);
        }
        AssertComplex(-(double)b / a, roots[0] + roots[1]);
        AssertComplex((double)c / a, roots[0] * roots[1]);
    }

    public static IEnumerable<object[]> Cubics()
    {
        float[] coefficients = { -6, -3, -1, 0, 1, 3, 6 };
        foreach (float a in coefficients) foreach (float b in coefficients) foreach (float c in coefficients)
            yield return new object[] { a, b, c };
        foreach (float scale in new[] { 1e-10f, 1e-3f, 1e3f, 1e10f })
        {
            yield return new object[] { 0f, 0f, -scale * scale * scale };
            yield return new object[] { -6 * scale, 11 * scale * scale, -6 * scale * scale * scale };
            yield return new object[] { scale, -1f, -scale };
        }
        foreach (float scale in new[] { 1e-30f, 1e-20f, 1e20f, 1e30f })
        {
            yield return new object[] { 0f, scale, 1f };
            yield return new object[] { 0f, 1f, scale };
            yield return new object[] { scale, scale, scale };
            yield return new object[] { -scale, scale, -scale };
        }
    }

    [Theory, MemberData(nameof(Cubics))]
    public void CubicRootsSatisfyAllVietaIdentitiesAndScaledResiduals(float a, float b, float c)
    {
        var roots = Maths.Cubic(a, b, c).Select(z => (Complex)z).ToArray();
        Assert.Equal(3, roots.Length);
        foreach (var root in roots)
        {
            double magnitude = root.Magnitude;
            double scale = ((magnitude + Math.Abs(a)) * magnitude + Math.Abs(b)) * magnitude + Math.Abs(c);
            Assert.True(Complex.Abs(((root + a) * root + b) * root + c) <= 5e-6 * scale + 1e-44,
                $"Coefficients ({a}, {b}, {c}); root {root} has an excessive residual.");
        }
        // Bound cancellation in the sum by the magnitudes of the individual rounded roots.
        Complex sum = roots[0] + roots[1] + roots[2];
        Assert.True(Complex.Abs(sum + a) <= 5e-6 * roots.Sum(z => z.Magnitude) + 1e-44);
        var pairs = new[] { roots[0] * roots[1], roots[0] * roots[2], roots[1] * roots[2] };
        double pairScale = pairs.Sum(z => z.Magnitude);
        Assert.True(Complex.Abs(pairs[0] + pairs[1] + pairs[2] - b) <= 5e-6 * pairScale + 1e-44);
        AssertComplex(-c, roots[0] * roots[1] * roots[2]);
    }

    [Fact]
    public void QuadraticRequiresANonzeroLeadingCoefficient() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Quadratic(0, 1, 1));

    private static void AssertComplex(double expected, Complex actual) =>
        Assert.True(Complex.Abs(actual - expected) <= 5e-6 * Math.Max(Math.Abs(expected), actual.Magnitude) + 1e-12,
            $"Expected {expected:G17}; actual {actual}.");

    private static void Rounded(double expected, float actual)
    {
        float rounded = (float)expected;
        if (float.IsNaN(rounded)) Assert.True(float.IsNaN(actual));
        else if (float.IsInfinity(rounded)) Assert.Equal(rounded, actual);
        else Close(expected, actual, 2.0 * float.Epsilon, 2e-6);
    }
}
