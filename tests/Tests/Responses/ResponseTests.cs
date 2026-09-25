using System.Numerics;
using UMapx.Core;
using UMapx.Response;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Response")]
public class ResponseTests
{
    static (IResponse filter, float[] b, float[] a) Create(string name)
    {
        if (name == "FIRCustom")
        {
            float[] b =
            {
                .25f,
                -.5f,
                .75f
            };
            return (new FIR(b), b, Array.Empty<float>());
        }

        if (name == "IIRCustom")
        {
            float[] b =
            {
                .25f,
                -.5f,
                .75f
            }, a =
            {
                .2f,
                .3f,
                -.1f
            };
            return (new IIR(b, a), b, a);
        }

        Type type = name.StartsWith("FIR") ? typeof(FIR) : typeof(IIR);
        var f = (IResponse)type.GetProperty(name[3..])!.GetValue(null)!;
        return (f, (float[])type.GetProperty("B")!.GetValue(f)!, type == typeof(IIR) ? (float[])type.GetProperty("A")!.GetValue(f)! : Array.Empty<float>());
    }

    public static IEnumerable<object[]> Cases()
    {
        foreach (string type in new[]
        {
            "FIR",
            "IIR"
        }

        )
            foreach (string name in new[]
            {
                "Custom",
                "LowPass",
                "HighPass",
                "BandPass",
                "Notch"
            }

            )
                foreach (bool complex in new[]
                {
                    false,
                    true
                }

                )
                    yield return new object[]
                    {
                        type + name,
                        complex
                    };
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void ReactionsSatisfyTheCausalDifferenceEquation(string name, bool complex)
    {
        var (f, b, a) = Create(name);
        var x = Enumerable.Range(0, 64).Select(i => new Complex32((float)Math.Sin(.19 * i), complex ? (float)Math.Cos(.13 * i) : 0)).ToArray();
        var expected = new Complex[x.Length];
        for (int n = 0; n < x.Length; n++)
        {
            for (int k = 0; k < b.Length && k <= n; k++)
                expected[n] += b[k] * (Complex)x[n - k];
            for (int k = 1; k < a.Length && k <= n; k++)
                expected[n] += a[k] * expected[n - k];
            expected[n] /= 1 - (a.Length > 0 ? a[0] : 0);
        }

        if (complex)
        {
            var actual = f.Reaction(x);
            for (int i = 0; i < x.Length; i++)
                Close(expected[i], actual[i], 2e-4);
        }
        else
        {
            var actual = f.Reaction(x.Select(z => z.Real).ToArray());
            for (int i = 0; i < x.Length; i++)
                Close(expected[i].Real, actual[i], 2e-4);
        }
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void FrequencyResponsesMatchTheTransferPolynomial(string name, bool complex)
    {
        var (f, b, a) = Create(name);
        var w = new[]
        {
            .17f,
            .51f,
            1.3f,
            2.7f
        };
        var z = w.Select(v => new Complex32(v, complex ? .03f : 0)).ToArray();
        var amplitudes = complex ? f.Amplitude(z) : f.Amplitude(w).Select(v => new Complex32(v, 0)).ToArray();
        var phases = complex ? f.Phase(z) : f.Phase(w).Select(v => new Complex32(v, 0)).ToArray();
        for (int i = 0; i < w.Length; i++)
        {
            Complex numerator = 0, denominator = 1;
            for (int k = 0; k < b.Length; k++)
                numerator += b[k] * Complex.Exp(-Complex.ImaginaryOne * (Complex)z[i] * k);
            for (int k = 0; k < a.Length; k++)
                denominator -= a[k] * Complex.Exp(-Complex.ImaginaryOne * (Complex)z[i] * k);
            double amplitude = (numerator / denominator).Magnitude, phase = numerator.Phase - denominator.Phase;
            Close(amplitude, amplitudes[i].Real);
            Close(phase, phases[i].Real);
            Close(0, amplitudes[i].Imag);
            Close(0, phases[i].Imag);
            if (complex)
            {
                Close(amplitude, f.Amplitude(z[i]).Real);
                Close(phase, f.Phase(z[i]).Real);
            }
            else
            {
                Close(amplitude, f.Amplitude(w[i]));
                Close(phase, f.Phase(w[i]));
            }
        }
    }

    [Theory]
    [InlineData(0f, .5f, true)]
    [InlineData(0f, 1.1f, false)]
    [InlineData(.9f, -.8f, false)]
    [InlineData(-.9f, 1.1f, true)]
    public void StabilityMatchesThePoleOfTheActualDifferenceEquation(float a0, float a1, bool stable)
    {
        Assert.Equal(stable, new IIR(new[] { 1f }, new[] { a0, a1 }).Stability);
    }

    [Fact]
    public void ReactionRejectsASingularInstantaneousFeedbackCoefficient()
    {
        var f = new IIR(new[] { 1f }, new[] { 1f });
        Assert.Throws<InvalidOperationException>(() => f.Reaction(new[] { 1f }));
        Assert.Throws<InvalidOperationException>(() => f.Reaction(new[] { Complex32.One }));
    }

    public static IEnumerable<object[]> PoleCases()
    {
        foreach (float a0 in new[]
        {
            -.9f,
            0f,
            .9f
        }

        )
            foreach (float pole in new[]
            {
                -1.05f,
                -.75f,
                0f,
                .5f,
                .98f,
                1f,
                1.05f
            }

            )
                foreach (bool paired in new[]
                {
                    false,
                    true
                }

                )
                    yield return new object[]
                    {
                        a0,
                        pole,
                        paired
                    };
    }

    [Theory, MemberData(nameof(PoleCases))]
    public void IirStabilityMatchesKnownRealAndComplexPoles(float a0, float radius, bool paired)
    {
        double[] polynomial = paired ? new[]
        {
            1.0,
            -(double)radius,
            (double)radius * radius
        }

        : new[]
        {
            1.0,
            -(double)radius
        };
        // The quadratic has conjugate poles radius * exp(+/- i*pi/3).
        var feedback = new float[polynomial.Length];
        feedback[0] = a0;
        for (int i = 1; i < feedback.Length; i++)
            feedback[i] = (float)(-(1.0 - a0) * polynomial[i]);
        var filter = new IIR
        {
            A = feedback,
            B = new[]
            {
                1f
            }
        };
        Assert.Equal(Math.Abs(radius) < 1, filter.Stability);
    }

    [Fact]
    public void IirHandlesRepeatedPolesZeroOrderAndUndefinedFeedback()
    {
        // (z - 1/2)^4, all roots strictly inside the unit disk.
        Assert.True(new IIR { A = new[] { 0f, 2f, -1.5f, .5f, -.0625f } }.Stability);
        Assert.True(new IIR { A = new[] { 2f } }.Stability);
        Assert.False(new IIR { A = new[] { 1f } }.Stability);
        Assert.False(new IIR { A = new[] { 0f, float.NaN } }.Stability);
        Assert.False(new IIR { A = new[] { float.PositiveInfinity } }.Stability);
    }

    [Fact]
    public void IirStabilityUsesTheSamePolynomialAsReaction() => Assert.False(new IIR(new[] { 1f }, new[] { .5f, -.75f }).Stability);
    [Fact]
    public void FirImpulseResponseEqualsCoefficients()
    {
        var fir = new FIR(new[] { 1f, 2f, 3f });
        Close(new[] { 1f, 2f, 3f, 0f }, fir.Reaction(new[] { 1f, 0f, 0f, 0f }));
    }
}
