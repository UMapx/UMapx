using System.Reflection;
using System.Text.Json;
using UMapx.Distribution;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Distribution")]
public class DistributionRepairTests
{
    public static IEnumerable<object[]> References()
    {
        using var stream = typeof(DistributionRepairTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.distribution-repair.json");
        using var document = JsonDocument.Parse(stream!);
        foreach (var value in document.RootElement.GetProperty("cases").EnumerateArray())
            yield return new object[] { value.GetProperty("name").GetString()!, value.GetProperty("member").GetString()!, value.GetRawText() };
    }

    [Theory]
    [MemberData(nameof(References))]
    public void RepairedMembersMatchIndependentHighPrecisionReferences(string name, string member, string json)
    {
        using var document = JsonDocument.Parse(json);
        var value = document.RootElement;
        var type = typeof(IDistribution).Assembly.GetType("UMapx.Distribution." + name)!;
        var signature = value.GetProperty("signature").EnumerateArray().Select(v => v.GetString() == "i" ? typeof(int) : typeof(float)).ToArray();
        var args = value.GetProperty("constructor").EnumerateArray().Select((v, i) => signature[i] == typeof(int) ? (object)(int)v.GetSingle() : v.GetSingle()).ToArray();
        var instance = type.GetConstructor(signature)!.Invoke(args);
        double actual;
        try
        {
            actual = Convert.ToDouble(member is "Function" or "Distribution"
                ? type.GetMethod(member, new[] { typeof(float) })!.Invoke(instance, new object[] { value.GetProperty("x").GetSingle() })
                : type.GetProperty(member)!.GetValue(instance));
        }
        catch (TargetInvocationException error) when (error.InnerException != null)
        {
            System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(error.InnerException).Throw();
            throw;
        }
        var reference = value.GetProperty("expected");
        double expected = reference.ValueKind == JsonValueKind.String
            ? double.Parse(reference.GetString()!, System.Globalization.CultureInfo.InvariantCulture)
            : reference.GetDouble();
        if (double.IsNaN(expected)) Assert.True(double.IsNaN(actual));
        else if (float.IsInfinity((float)expected)) Assert.Equal((double)(float)expected, actual);
        else if (member == "Median" && name is "Poisson" or "Binomial") Assert.Equal(expected, actual);
        // Retain relative accuracy in small probabilities instead of accepting zero under a large absolute budget.
        else Close(expected, actual, 2 * (double)float.Epsilon, 2e-5);
    }

    public static IEnumerable<object[]> ModeCases()
    {
        foreach (float scale in new[] { .1f, 1f, 7f })
            foreach (float ratio in new[] { 0f, .5f, 1f, 1.001f, 1.1f, 2f, 10f })
                foreach (int sign in new[] { -1, 1 }) yield return new object[] { scale, ratio, sign };
    }

    [Theory]
    [MemberData(nameof(ModeCases))]
    public void FoldedNormalModeSolvesTheDensityStationarityEquation(float scale, float ratio, int sign)
    {
        float mu = sign * ratio * scale;
        var law = new FoldedNormal(mu, scale);
        double mode = Assert.Single(law.Mode), a = Math.Abs((double)mu), s = scale;
        Assert.InRange(mode, 0, a);
        if (a <= s) Assert.Equal(0, mode);
        else
        {
            Assert.True(mode > 0);
            Close(mode, a * Math.Tanh(a * mode / (s * s)), 1e-7 * s, 2e-6);
        }
        double Density(double x) => (Math.Exp(-.5 * Math.Pow((x - a) / s, 2)) + Math.Exp(-.5 * Math.Pow((x + a) / s, 2))) / (s * Math.Sqrt(2 * Math.PI));
        foreach (double x in new[] { 0, Math.Max(0, mode - .1 * s), mode + .1 * s, a, a + s })
            Assert.True(Density(mode) + 1e-12 / s >= Density(x));
    }

    [Theory]
    [InlineData(.0001f)] [InlineData(.1f)] [InlineData(.7f)] [InlineData(1f)] [InlineData(3f)] [InlineData(100f)]
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
    [InlineData(.001f)] [InlineData(.1f)] [InlineData(1f)] [InlineData(20f)] [InlineData(10000f)]
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
    public void RepairedProbabilityLawsHaveDefinedEndpointsAndPropagateNaN()
    {
        IDistribution[] laws = { new Poisson(.9f), new Binomial(5, .3f), new PowerNormal(.1f), new PowerLognormal(.1f, .7f), new FisherZ(4, 12),
            new Gaussian(1.25f, .5f), new GeneralizedNormal(.5f, 1.5f, 2.5f) };
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
}
