using System.Numerics;
using System.Text.Json;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class ArithmeticTests
{
    private static readonly Dictionary<string, Func<double, double>> RealFunctions = new()
    {
        ["Pow"] = x => x * x,
        ["Exp"] = Math.Exp,
        ["Log"] = Math.Log,
        ["Log10"] = Math.Log10,
        ["Log2"] = Math.Log2,
        ["Sqrt"] = Math.Sqrt,
        ["Abs"] = Math.Abs,
        ["Floor"] = Math.Floor,
        ["Ceil"] = Math.Ceiling,
        ["Round"] = Math.Round,
        ["Cos"] = Math.Cos,
        ["Sin"] = Math.Sin,
        ["Tan"] = Math.Tan,
        ["Ctan"] = x => 1 / Math.Tan(x),
        ["Sec"] = x => 1 / Math.Cos(x),
        ["Cosc"] = x => 1 / Math.Sin(x),
        ["Asin"] = Math.Asin,
        ["Acos"] = Math.Acos,
        ["Atan"] = Math.Atan,
        ["Actan"] = x => Math.PI / 2 - Math.Atan(x),
        ["Asec"] = x => Math.Acos(1 / x),
        ["Acosc"] = x => Math.Asin(1 / x),
        ["Sinh"] = Math.Sinh,
        ["Cosh"] = Math.Cosh,
        ["Tanh"] = Math.Tanh,
        ["Ctanh"] = x => 1 / Math.Tanh(x),
        ["Sech"] = x => 1 / Math.Cosh(x),
        ["Cosch"] = x => 1 / Math.Sinh(x),
        ["Asinh"] = Math.Asinh,
        ["Acosh"] = Math.Acosh,
        ["Atanh"] = Math.Atanh,
        ["Actanh"] = x => Math.Atanh(1 / x),
        ["Asech"] = x => Math.Acosh(1 / x),
        ["Acosch"] = x => Math.Asinh(1 / x)
    };
    public static IEnumerable<object[]> RealCases()
    {
        foreach (var pair in RealFunctions)
            foreach (float x in new[]
            {
                -10000f,
                -100,
                -5,
                -2,
                -.9f,
                -.1f,
                .1f,
                .5f,
                .9f,
                1,
                2,
                5,
                100,
                10000
            }

            )
            {
                double expected = pair.Value(x);
                if (double.IsFinite(expected) && Math.Abs(expected) <= float.MaxValue)
                    yield return new object[]
                    {
                        pair.Key,
                        x,
                        expected
                    };
            }
    }

    [Theory, MemberData(nameof(RealCases))]
    public void RealElementaryFunctionsAgreeWithDoublePrecision(string name, float x, double expected)
    {
        var method = typeof(Maths).GetMethod(name, new[] { typeof(float) })!;
        Close(expected, Convert.ToDouble(method.Invoke(null, new object[] { x })), 3e-6, 3e-5);
    }

    private static Complex Asinh(Complex z) => z.Real < 0 ? -Asinh(-z) : Complex.Log(z + Complex.Sqrt(z * z + 1));
    private static Complex Acosh(Complex z) => Complex.Log(z + Complex.Sqrt(z - 1) * Complex.Sqrt(z + 1));
    private static readonly Dictionary<string, Func<Complex, Complex>> ComplexFunctions = new()
    {
        ["Exp"] = Complex.Exp,
        ["Log"] = Complex.Log,
        ["Log10"] = Complex.Log10,
        ["Log2"] = z => Complex.Log(z) / Math.Log(2),
        ["Sqrt"] = Complex.Sqrt,
        ["Sin"] = Complex.Sin,
        ["Cos"] = Complex.Cos,
        ["Tan"] = Complex.Tan,
        ["Ctan"] = z => 1 / Complex.Tan(z),
        ["Sec"] = z => 1 / Complex.Cos(z),
        ["Cosc"] = z => 1 / Complex.Sin(z),
        ["Asin"] = Complex.Asin,
        ["Acos"] = Complex.Acos,
        ["Atan"] = Complex.Atan,
        ["Actan"] = z => Complex.Atan(1 / z),
        ["Asec"] = z => Complex.Acos(1 / z),
        ["Acosc"] = z => Complex.Asin(1 / z),
        ["Sinh"] = Complex.Sinh,
        ["Cosh"] = Complex.Cosh,
        ["Tanh"] = Complex.Tanh,
        ["Ctanh"] = z => 1 / Complex.Tanh(z),
        ["Sech"] = z => 1 / Complex.Cosh(z),
        ["Cosch"] = z => 1 / Complex.Sinh(z),
        ["Asinh"] = Asinh,
        ["Acosh"] = Acosh,
        ["Atanh"] = z => (Complex.Log(1 + z) - Complex.Log(1 - z)) / 2,
        ["Actanh"] = z => (Complex.Log(1 + 1 / z) - Complex.Log(1 - 1 / z)) / 2,
        ["Asech"] = z => Acosh(1 / z),
        ["Acosch"] = z => Asinh(1 / z)
    };
    public static IEnumerable<object[]> ComplexCases()
    {
        foreach (string name in ComplexFunctions.Keys)
            foreach (var z in new[]
            {
                new Complex32(.2f, .4f),
                new Complex32(-2, .5f),
                new Complex32(-2, -.5f),
                new Complex32(4, 1),
                new Complex32(1, 20),
                new Complex32(-1, -20)
            }

            )
                yield return new object[]
                {
                    name,
                    z.Real,
                    z.Imag
                };
    }

    [Theory, MemberData(nameof(ComplexCases))]
    public void ComplexElementaryFunctionsRespectPrincipalValues(string name, float real, float imaginary)
    {
        var z = new Complex32(real, imaginary);
        var method = typeof(Maths).GetMethod(name, new[] { typeof(Complex32) })!;
        Close(ComplexFunctions[name](z), (Complex32)method.Invoke(null, new object[] { z })!, 5e-5, 8e-5);
    }

    [Theory]
    [InlineData(7)]
    [InlineData(53)]
    [InlineData(731)]
    public void ComplexArithmeticAndConversionsAgreeWithSystemNumerics(int seed)
    {
        var random = new Random(seed);
        for (int i = 0; i < 40; i++)
        {
            var a = new Complex32((float)(random.NextDouble() * 4 - 2), (float)(random.NextDouble() * 4 - 2));
            var b = new Complex32((float)(random.NextDouble() + .1), (float)(random.NextDouble() + .1));
            float s = .25f + (float)random.NextDouble();
            Complex x = a, y = b;
            Close(x + y, a + b);
            Close(x - y, a - b);
            Close(x * y, a * b);
            Close(x / y, a / b);
            Close(x + s, a + s);
            Close(s + x, s + a);
            Close(x - s, a - s);
            Close(s - x, s - a);
            Close(x * s, a * s);
            Close(s * x, s * a);
            Close(x / s, a / s);
            Close(s / x, s / a);
            Close(-x, -a);
            Close(x, +a);
            Close(Complex.Conjugate(x), a.Conjugate);
            Close(x.Magnitude, a.Abs);
            Close(x.Magnitude * x.Magnitude, a.Abs2);
            Close(x.Phase, Maths.Angle(a));
            Close(x.Magnitude, Maths.Abs(a));
            Close(x, Maths.FromPolar((float)x.Magnitude, (float)x.Phase));
            Close(x, Complex32.FromPolarCoordinates((float)x.Magnitude, (float)x.Phase));
            Close(Complex.Pow(x, s), Maths.Pow(a, s));
            Close(Complex.Pow(x, y), Maths.Pow(a, b));
            Close(Complex.Pow(s, y), Maths.Pow(s, b));
            Close(Complex.Pow(x, 1 / s), Maths.Sqrt(a, s));
            Close(Complex.Pow(x, 1 / y), Maths.Sqrt(a, b));
            Close(Complex.Log(x) / Math.Log(s), Maths.Log(a, s));
            Close(new Complex(Math.Round(a.Real), Math.Round(a.Imag)), Maths.Round(a));
            Close(new Complex(Math.Round(a.Real, 2), Math.Round(a.Imag, 2)), Maths.Round(a, 2));
            Assert.Equal(a, a.Clone());
            Assert.True(a == a.Clone());
            Assert.False(a != a.Clone());
            Assert.False(Complex32.IsNaN(a));
            Assert.False(Complex32.IsInfinity(a));
            Assert.Equal(a.GetHashCode(), a.Clone().GetHashCode());
        }

        Assert.True(Complex32.IsNaN(new Complex32(float.NaN, 0)));
        Assert.True(Complex32.IsInfinity(new Complex32(0, float.PositiveInfinity)));
    }

    private static Quaternion Q(Quaternion32 value) => new(value.X, value.Y, value.Z, value.W);
    private static Quaternion32 Q(Quaternion value) => new(value.X, value.Y, value.Z, value.W);
    private static void EqualQuaternion(Quaternion expected, Quaternion32 actual)
    {
        Close(expected.X, actual.X);
        Close(expected.Y, actual.Y);
        Close(expected.Z, actual.Z);
        Close(expected.W, actual.W);
    }

    [Theory]
    [InlineData(0f)]
    [InlineData(.2f)]
    [InlineData(.75f)]
    [InlineData(1f)]
    public void QuaternionAlgebraAndInterpolationAgreeWithSystemNumerics(float amount)
    {
        var a = Quaternion32.FromYPR(.2f, -.4f, .7f);
        var b = Quaternion32.FromYPR(-.6f, .3f, -.8f);
        EqualQuaternion(Quaternion.CreateFromYawPitchRoll(.2f, -.4f, .7f), a);
        EqualQuaternion(Q(a) + Q(b), a + b);
        EqualQuaternion(Q(a) - Q(b), a - b);
        EqualQuaternion(Q(a) * Q(b), a * b);
        EqualQuaternion(Q(a) / Q(b), a / b);
        EqualQuaternion(-Q(a), -a);
        EqualQuaternion(Q(a) * 2, a * 2);
        EqualQuaternion(Q(a) * .5f, a / 2);
        EqualQuaternion(Quaternion.Conjugate(Q(a)), a.Conjugate);
        EqualQuaternion(Quaternion.Inverse(Q(a)), a.Inverse);
        EqualQuaternion(Quaternion.Normalize(Q(a) * 2), (a * 2).Normalize);
        EqualQuaternion(Quaternion.Concatenate(Q(a), Q(b)), Quaternion32.Concatenate(a, b));
        EqualQuaternion(Quaternion.Lerp(Q(a), Q(b), amount), Quaternion32.Lerp(a, b, amount));
        EqualQuaternion(Quaternion.Slerp(Q(a), Q(b), amount), Quaternion32.Slerp(a, b, amount));
        Close(Quaternion.Dot(Q(a), Q(b)), Quaternion32.Dot(a, b));
        Close(Q(a).Length(), a.Abs);
        Close(Q(a).LengthSquared(), a.SquaredAbs);
        Assert.True(Quaternion32.Identity.IsIdentity);
        Assert.Equal(a, a.Clone());
        Assert.True(a == a.Clone());
        Assert.False(a != a.Clone());
        Assert.Equal(a.GetHashCode(), a.Clone().GetHashCode());
    }

    [Fact]
    public void ScalarRangesAndRoundingFollowTheirDefinitions()
    {
        foreach (float x in new[]
        {
            -300f,
            -1.25f,
            0,
            1.25f,
            100,
            300
        }

        )
        {
            Close(Math.Clamp(x, -2, 2), Maths.Range(x, -2, 2));
            Assert.Equal(x >= -2 && x <= 2, Maths.IsRange(x, -2, 2));
            Close((x + 2) / 4, Maths.Normalize(x, -2, 2));
            Assert.Equal((byte)Math.Clamp((int)x, 0, 255), Maths.Byte(x));
            Assert.Equal((sbyte)Math.Clamp((int)x, -128, 127), Maths.sByte(x));
            Close(Math.Round(x, 1), Maths.Round(x, 1));
        }

        foreach (int x in new[]
        {
            -300,
            -1,
            0,
            1,
            100,
            300
        }

        )
        {
            Assert.Equal(Math.Clamp(x, -2, 2), Maths.Range(x, -2, 2));
            Assert.Equal(x >= -2 && x <= 2, Maths.IsRange(x, -2, 2));
            Assert.Equal((byte)Math.Clamp(x, 0, 255), Maths.Byte(x));
            Assert.Equal((sbyte)Math.Clamp(x, -128, 127), Maths.sByte(x));
        }

        foreach (float x in new[]
        {
            .25f,
            1,
            4,
            16
        }

        )
        {
            Close(Math.Log(x, 3), Maths.Log(x, 3));
            Close(Math.Pow(x, 1.0 / 3), Maths.Sqrt(x, 3));
            Close(x * x, Maths.Pow(x));
        }

        foreach (float a in new[]
        {
            -2f,
            .5f,
            3
        }

        )
            foreach (float b in new[]
            {
                -3f,
                .25f,
                2
            }

            )
            {
                Close(Math.Max(a, b), Maths.Max(a, b));
                Close(Math.Min(a, b), Maths.Min(a, b));
                Close(Math.Max(a, Math.Max(b, 1)), Maths.Max(a, b, 1));
                Close(Math.Min(a, Math.Min(b, 1)), Maths.Min(a, b, 1));
                Close(Math.Atan2(a, b), Maths.Atan2(a, b));
                Close(Math.Sqrt((double)a * a + (double)b * b), Maths.Hypotenuse(a, b));
            }
    }

    public static IEnumerable<object[]> References()
    {
        using var stream = typeof(ArithmeticTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.arithmetic.json")!;
        using var data = JsonDocument.Parse(stream);
        foreach (var row in data.RootElement.GetProperty("cases").EnumerateArray())
            yield return new object[]
            {
                row[0].GetString()!,
                row[1].GetSingle(),
                row[2].GetSingle(),
                row[3].GetSingle(),
                row[4].GetSingle(),
                row[5].GetDouble(),
                row[6].GetDouble()
            };
    }

    [Theory, MemberData(nameof(References))]
    public void ComplexBranchesAndPowersMatchIndependentHighPrecisionReferences(string operation, float re, float im, float exponentRe, float exponentIm, double expectedRe, double expectedIm)
    {
        var z = new Complex32(re, im);
        Complex32 actual = operation switch
        {
            "Asinh" => Maths.Asinh(z),
            "Acosh" => Maths.Acosh(z),
            "Actan" => Maths.Actan(z),
            "Tanh" => Maths.Tanh(z),
            "Ctanh" => Maths.Ctanh(z),
            "Sech" => Maths.Sech(z),
            "Cosch" => Maths.Cosch(z),
            "Pow" => Maths.Pow(re, new Complex32(exponentRe, exponentIm)),
            _ => throw new ArgumentException(operation)
        };
        Close(new Complex(expectedRe, expectedIm), actual, 4.0 * float.Epsilon, 3e-6);
    }

    public static IEnumerable<object[]> RealBoundaries()
    {
        foreach (float value in new[]
        {
            float.Epsilon,
            1e-30f,
            1e-10f,
            1e-4f,
            .5f,
            1f,
            1.0000001f,
            20f,
            89f,
            100f,
            1e20f,
            float.MaxValue
        }

        )
        {
            yield return new object[]
            {
                value
            };
            yield return new object[]
            {
                -value
            };
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
        foreach (float x in new[]
        {
            0f,
            negativeZero,
            float.PositiveInfinity,
            float.NegativeInfinity,
            float.NaN
        }

        )
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
        foreach (float scale in new[]
        {
            1e-30f,
            1e-10f,
            1f,
            1e10f,
            1e30f
        }

        )
            foreach (var pair in new[]
            {
                (0f, 0f),
                (1f, 1f),
                (0f, -1f),
                (0f, 1f),
                (-3f, 2f),
                (1e6f, 1f),
                (-1e6f, 1f)
            }

            )
                yield return new object[]
                {
                    scale,
                    pair.Item1 * scale,
                    pair.Item2 * scale
                };
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
        float[] coefficients =
        {
            -6,
            -3,
            -1,
            0,
            1,
            3,
            6
        };
        foreach (float a in coefficients)
            foreach (float b in coefficients)
                foreach (float c in coefficients)
                    yield return new object[]
                    {
                        a,
                        b,
                        c
                    };
        foreach (float scale in new[]
        {
            1e-10f,
            1e-3f,
            1e3f,
            1e10f
        }

        )
        {
            yield return new object[]
            {
                0f,
                0f,
                -scale * scale * scale
            };
            yield return new object[]
            {
                -6 * scale,
                11 * scale * scale,
                -6 * scale * scale * scale
            };
            yield return new object[]
            {
                scale,
                -1f,
                -scale
            };
        }

        foreach (float scale in new[]
        {
            1e-30f,
            1e-20f,
            1e20f,
            1e30f
        }

        )
        {
            yield return new object[]
            {
                0f,
                scale,
                1f
            };
            yield return new object[]
            {
                0f,
                1f,
                scale
            };
            yield return new object[]
            {
                scale,
                scale,
                scale
            };
            yield return new object[]
            {
                -scale,
                scale,
                -scale
            };
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
            Assert.True(Complex.Abs(((root + a) * root + b) * root + c) <= 5e-6 * scale + 1e-44, $"Coefficients ({a}, {b}, {c}); root {root} has an excessive residual.");
        }

        // Bound cancellation in the sum by the magnitudes of the individual rounded roots.
        Complex sum = roots[0] + roots[1] + roots[2];
        Assert.True(Complex.Abs(sum + a) <= 5e-6 * roots.Sum(z => z.Magnitude) + 1e-44);
        var pairs = new[]
        {
            roots[0] * roots[1],
            roots[0] * roots[2],
            roots[1] * roots[2]
        };
        double pairScale = pairs.Sum(z => z.Magnitude);
        Assert.True(Complex.Abs(pairs[0] + pairs[1] + pairs[2] - b) <= 5e-6 * pairScale + 1e-44);
        AssertComplex(-c, roots[0] * roots[1] * roots[2]);
    }

    [Fact]
    public void QuadraticRequiresANonzeroLeadingCoefficient() => Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Quadratic(0, 1, 1));
    private static void AssertComplex(double expected, Complex actual) => Assert.True(Complex.Abs(actual - expected) <= 5e-6 * Math.Max(Math.Abs(expected), actual.Magnitude) + 1e-12, $"Expected {expected:G17}; actual {actual}.");
    private static void Rounded(double expected, float actual)
    {
        float rounded = (float)expected;
        if (float.IsNaN(rounded))
            Assert.True(float.IsNaN(actual));
        else if (float.IsInfinity(rounded))
            Assert.Equal(rounded, actual);
        else
            Close(expected, actual, 2.0 * float.Epsilon, 2e-6);
    }

    [Theory]
    [InlineData(-1f)]
    [InlineData(0f)]
    [InlineData(1f)]
    public void QuadraticHandlesZeroLinearCoefficient(float constant)
    {
        var roots = Maths.Quadratic(1, 0, constant);
        Assert.Equal(2, roots.Length);
        foreach (var root in roots)
            Close(0, (Complex)root * (Complex)root + constant);
        Close(0, (Complex)roots[0] + (Complex)roots[1]);
        Close(constant, (Complex)roots[0] * (Complex)roots[1]);
    }

    [Theory]
    [InlineData(0f, 0f, -1f)]
    [InlineData(0f, -3f, -2f)]
    public void CubicUsesRealCubeRootsForNegativeRadicands(float a, float b, float c)
    {
        var roots = Maths.Cubic(a, b, c);
        Assert.Equal(3, roots.Length);
        foreach (Complex root in roots)
            Close(0, ((root + a) * root + b) * root + c);
    }

    [Fact]
    public void ComplexArccotangentAgreesWithPositiveRealBranch() => Close(new Complex(Math.PI / 4, 0), Maths.Actan(new Complex32(1, 0)));
    [Fact]
    public void ComplexAcoshUsesPrincipalBranchOnNegativeRealAxis() => Close(new Complex(Math.Acosh(2), Math.PI), Maths.Acosh(new Complex32(-2, 0)));
    [Fact]
    public void TanhAvoidsInfinityDividedByInfinity() => Close(1, Maths.Tanh(100f));
    [Fact]
    public void AsinhAvoidsCancellationOnNegativeArguments() => Close(Math.Asinh(-10000), Maths.Asinh(-10000f));
    [Theory]
    [InlineData(1e-20f)]
    [InlineData(1e20f)]
    public void ComplexDivisionIsInvariantToScaling(float scale)
    {
        var z = new Complex32(scale, scale);
        Close(Complex.One, z / z);
    }

    [Fact]
    public void NegativeRealBaseSupportsComplexExponent() => Close(Complex.ImaginaryOne, Maths.Pow(-1f, new Complex32(.5f, 0)));
}
