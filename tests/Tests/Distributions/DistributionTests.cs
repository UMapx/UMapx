using Poisson = UMapx.Distribution.Poisson;
using System.Text.Json;
using UMapx.Distribution;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Distribution")]
public class DistributionTests
{
    public static IEnumerable<object[]> Cases()
    {
        using var stream = typeof(DistributionTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.distribution-catalog.json");
        using var document = JsonDocument.Parse(stream!);
        foreach (var law in document.RootElement.EnumerateArray())
            foreach (string member in new[]
            {
                "Mode",
                "Median"
            }

            )
                yield return new object[]
                {
                    law.GetProperty("name").GetString()!,
                    member,
                    law.GetRawText()
                };
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void ModesAndMediansSatisfyTheirProbabilityDefinitions(string name, string member, string json)
    {
        using var document = JsonDocument.Parse(json);
        var law = document.RootElement;
        var type = typeof(IDistribution).Assembly.GetType("UMapx.Distribution." + name)!;
        var signature = law.GetProperty("signature").EnumerateArray().Select(v => v.GetString() == "i" ? typeof(int) : typeof(float)).ToArray();
        var args = law.GetProperty("constructor").EnumerateArray().Select((v, i) => signature[i] == typeof(int) ? (object)(int)v.GetSingle() : v.GetSingle()).ToArray();
        var d = (IDistribution)type.GetConstructor(signature)!.Invoke(args);
        if (law.GetProperty("unsupported").EnumerateArray().Any(v => v.GetString() == member))
        {
            // An explicit unsupported getter is an API limitation, not a successful numerical calculation.
            var exception = Record.Exception(() => type.GetProperty(member)!.GetValue(d));
            if (exception != null)
                Assert.IsType<NotSupportedException>(exception.InnerException ?? exception);
            else
                Assert.True(float.IsNaN(Convert.ToSingle(type.GetProperty(member)!.GetValue(d))));
            return;
        }

        double Evaluate(string method, double x) => Convert.ToDouble(type.GetMethod(method, new[] { typeof(float) })!.Invoke(d, new object[] { (float)x }));
        bool discrete = law.GetProperty("discrete").GetBoolean();
        if (member == "Median")
        {
            float median = d.Median;
            Assert.True(float.IsFinite(median), $"Expected a finite median; actual {median}.");
            if (law.GetProperty("unsupported").EnumerateArray().Any(v => v.GetString() == "Distribution"))
                return;
            double right = Evaluate("Distribution", median), left = discrete ? Evaluate("Distribution", Math.Ceiling(median) - 1) : right;
            Assert.True(left <= .501 && right >= .499, $"Median {median:G9}: F(m-)={left:G9}, F(m)={right:G9}.");
        }
        else
        {
            float[] modes = d.Mode;
            Assert.NotEmpty(modes);
            foreach (float mode in modes)
            {
                if (float.IsNaN(mode))
                {
                    // These APIs explicitly represent a flat mode set with NaN.
                    Assert.Contains(name, new[] { "Uniform", "UniformDiscrete", "Trapezoidal" });
                    continue;
                }

                Assert.InRange(mode, d.Support.Min, d.Support.Max);
                if (!discrete && (mode == d.Support.Min || mode == d.Support.Max))
                    continue; // Endpoint density conventions vary.
                double center = Evaluate("Function", mode);
                Assert.False(double.IsNaN(center));
                double step = discrete ? 1 : .02 * (1 + Math.Abs(mode));
                foreach (double x in new[]
                {
                    mode - step,
                    mode + step
                }

                )
                    if (x > d.Support.Min && x < d.Support.Max)
                    {
                        double neighbor = Evaluate("Function", x);
                        Assert.True(center + 2e-5 + Math.Abs(center) * 2e-5 >= neighbor, $"Mode {mode:G9}: density {center:G9} < density {neighbor:G9} at {x:G9}.");
                    }
            }
        }
    }

    [Theory]
    [InlineData(0f, 1f, 1f)]
    [InlineData(.5f, 2f, .7f)]
    public void BirnbaumSaundersModeIsThePositiveStationaryPoint(float location, float scale, float shape)
    {
        // Differentiate log f(x): t^3+(1+g^2)t^2+(3g^2-1)t-1=0, t=(x-location)/scale.
        double lo = 0, hi = 1, g2 = (double)shape * shape;
        for (int i = 0; i < 80; i++)
        {
            double t = (lo + hi) / 2;
            double value = t * t * t + (1 + g2) * t * t + (3 * g2 - 1) * t - 1;
            if (value > 0)
                hi = t;
            else
                lo = t;
        }

        var mode = new BirnbaumSaunders(location, scale, shape).Mode;
        Assert.Single(mode);
        Close(location + scale * (lo + hi) / 2, mode[0], 3e-5);
    }

    [Fact]
    public void FoldedNormalWithLocationSmallerThanScaleHasModeZero()
    {
        // A positive stationary point would satisfy x=mu*tanh(mu*x/sigma^2).
        // tanh(t)<t excludes such a point when |mu|<=sigma.
        var modes = new FoldedNormal(.75f, 1.25f).Mode;
        Assert.Single(modes);
        Close(0, modes[0]);
    }

    public static IEnumerable<object[]> ModeCases()
    {
        foreach (float scale in new[]
        {
            .1f,
            1f,
            7f
        }

        )
            foreach (float ratio in new[]
            {
                0f,
                .5f,
                1f,
                1.001f,
                1.1f,
                2f,
                10f
            }

            )
                foreach (int sign in new[]
                {
                    -1,
                    1
                }

                )
                    yield return new object[]
                    {
                        scale,
                        ratio,
                        sign
                    };
    }

    [Theory]
    [MemberData(nameof(ModeCases))]
    public void FoldedNormalModeSolvesTheDensityStationarityEquation(float scale, float ratio, int sign)
    {
        float mu = sign * ratio * scale;
        var law = new FoldedNormal(mu, scale);
        double mode = Assert.Single(law.Mode), a = Math.Abs((double)mu), s = scale;
        Assert.InRange(mode, 0, a);
        if (a <= s)
            Assert.Equal(0, mode);
        else
        {
            Assert.True(mode > 0);
            Close(mode, a * Math.Tanh(a * mode / (s * s)), 1e-7 * s, 2e-6);
        }

        double Density(double x) => (Math.Exp(-.5 * Math.Pow((x - a) / s, 2)) + Math.Exp(-.5 * Math.Pow((x + a) / s, 2))) / (s * Math.Sqrt(2 * Math.PI));
        foreach (double x in new[]
        {
            0,
            Math.Max(0, mode - .1 * s),
            mode + .1 * s,
            a,
            a + s
        }

        )
            Assert.True(Density(mode) + 1e-12 / s >= Density(x));
    }

    [Theory]
    [InlineData(.0001f)]
    [InlineData(.1f)]
    [InlineData(.7f)]
    [InlineData(1f)]
    [InlineData(3f)]
    [InlineData(100f)]
    public void BirnbaumSaundersModeIsThePositiveStationaryMaximum(float shape)
    {
        var law = new BirnbaumSaunders(0, 1, shape);
        double t = Assert.Single(law.Mode), g = shape;
        Assert.InRange(t, double.Epsilon, 1);
        Close(0, t * t * t + (1 + g * g) * t * t + (3 * g * g - 1) * t - 1, 4e-7, 0);
        double LogDensity(double u) => Math.Log(u + 1) - 1.5 * Math.Log(u) - (u + 1 / u - 2) / (2 * g * g);
        Assert.True(LogDensity(t) >= LogDensity(t * .99));
        Assert.True(LogDensity(t) >= LogDensity(t * 1.01));
    }

    [Theory]
    [InlineData(.001f)]
    [InlineData(.1f)]
    [InlineData(1f)]
    [InlineData(20f)]
    [InlineData(10000f)]
    public void GompertzMomentsRespectRateScaling(float shape)
    {
        var unit = new Gompertz(shape, 1);
        var scaled = new Gompertz(shape, 4);
        Assert.True(unit.Mean > 0 && unit.Variance > 0);
        Close(unit.Mean / 4, scaled.Mean);
        Close(unit.Variance / 16, scaled.Variance);
    }

    [Fact]
    public void DiscretePointMassAndModeDoNotRoundOntoDifferentIntegers()
    {
        // 16777217 is not representable as binary32; the adjacent integer is outside this point mass.
        var degenerate = new Binomial(16777217, 1);
        Assert.Equal(0, degenerate.Function(16777216f));
        Assert.Equal(0, degenerate.Distribution(16777216f));
        Assert.Equal(0, degenerate.Function(16777218f));
        Assert.Equal(1, degenerate.Distribution(16777218f));
        Assert.Equal(new[] { 0f, 1f }, new Poisson(1).Mode);
        Assert.Equal(new[] { 1f }, new Poisson(MathF.BitIncrement(1)).Mode);
        Assert.Equal(new[] { 0f }, new Poisson(MathF.BitDecrement(1)).Mode);
    }

    [Fact]
    public void ProbabilityLawsHaveDefinedEndpointsAndPropagateNaN()
    {
        IDistribution[] laws =
        {
            new Poisson(.9f),
            new Binomial(5, .3f),
            new PowerNormal(.1f),
            new PowerLognormal(.1f, .7f),
            new FisherZ(4, 12),
            new Gaussian(1.25f, .5f),
            new GeneralizedNormal(.5f, 1.5f, 2.5f)
        };
        foreach (dynamic law in laws)
        {
            Assert.Equal(0, law.Distribution(float.NegativeInfinity));
            Assert.Equal(1, law.Distribution(float.PositiveInfinity));
            Assert.Equal(0, law.Function(float.NegativeInfinity));
            Assert.Equal(0, law.Function(float.PositiveInfinity));
            Assert.True(float.IsNaN(law.Function(float.NaN)));
            Assert.True(float.IsNaN(law.Distribution(float.NaN)));
        }

        Assert.Equal(0, new PowerLognormal(.1f, .7f).Function(0));
        Assert.Equal(0, new PowerLognormal(.1f, .7f).Distribution(0));
        var tukey = new TukeyLambda();
        Assert.Equal(0, tukey.Function(float.NegativeInfinity));
        Assert.Equal(0, tukey.Function(float.PositiveInfinity));
        Assert.True(float.IsNaN(tukey.Function(float.NaN)));
    }

    [Theory]
    [InlineData(-2f)]
    [InlineData(-.5f)]
    [InlineData(.5f)]
    [InlineData(2f)]
    public void TimeFrequencyKernelsMatchTheirFourierPair(float tau)
    {
        var choi = new ChoiWilliams(.3f);
        var cone = new ConeShape(.3f);
        foreach (float eta in new[]
        {
            -2f,
            -.1f,
            0f,
            .5f,
            2f
        }

        )
        {
            Close(Math.Exp(-.3f * eta * eta * tau * tau), choi.Function(eta, tau));
            double v = Math.PI * eta * tau;
            Close((v == 0 ? 1 : Math.Sin(v) / v) * Math.Exp(-2 * Math.PI * .3f * tau * tau), cone.Function(eta, tau));
        }

        // Inverse Fourier transform of sinc(pi*eta*tau): an even rectangular pulse of width |tau|.
        foreach (float t in new[]
        {
            0f,
            .1f,
            3f
        }

        )
            Close(Math.Abs(t) < Math.Abs(tau) / 2 ? Math.Exp(-2 * Math.PI * .3f * tau * tau) / Math.Abs(tau) : 0, cone.Distribution(t, tau));
    }

    [Theory]
    [InlineData(.125f)]
    [InlineData(.5f)]
    [InlineData(2f)]
    public void ConeTimeKernelHasTheMassOfItsFourierValueAtZero(float tau)
    {
        var kernel = new ConeShape(.2f);
        foreach (float signed in new[]
        {
            -tau,
            tau
        }

        )
        {
            double sum = 0;
            for (int i = 0; i < 1000; i++)
                sum += kernel.Distribution((float)(-tau / 2.0 + (i + .5) * tau / 1000), signed) * tau / 1000;
            Close(kernel.Function(0, signed), sum, 2e-6, 2e-6);
            Close(0, kernel.Distribution(tau, signed), 0, 0);
        }
    }

    [Theory]
    [InlineData(0)]
    [InlineData(5)]
    public void BinomialCertainSuccessHasUnitMass(int n) => Close(1, new Binomial(n, 1).Function(n));
    [Fact]
    public void PoissonMassIsFiniteAtItsModeForLambda100() => Close(.0398609968091471, new Poisson(100).Function(100));
    [Fact]
    public void BinomialMedianSatisfiesBothHalfProbabilityInequalities() => Close(1, new Binomial(2, .7f).Median);
    [Fact]
    public void PoissonMedianCanBeOneBelowLambdaOne() => Close(1, new Poisson(.9f).Median);
    [Fact]
    public void InverseChiSquareEntropyAgreesWithInverseGamma() => Close(1 - Math.Log(2) + 2 * .57721566490153286, new InverseChiSquare(2).Entropy);
}
