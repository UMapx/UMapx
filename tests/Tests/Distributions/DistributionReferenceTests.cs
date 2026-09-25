using System.Reflection;
using System.Runtime.ExceptionServices;
using System.Text.Json;
using UMapx.Distribution;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Distribution")]
public class DistributionReferenceTests
{
    public static IEnumerable<object[]> Cases() => ReadCases("distributions.json");
    public static IEnumerable<object[]> EdgeCases() => ReadCases("distribution-edge-cases.json", true);
    private static IEnumerable<object[]> ReadCases(string resource, bool wrapped = false)
    {
        using var stream = typeof(DistributionReferenceTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data." + resource);
        using var document = JsonDocument.Parse(stream!);
        var rows = wrapped ? document.RootElement.GetProperty("cases") : document.RootElement;
        foreach (var value in rows.EnumerateArray())
            yield return new object[]
            {
                value.GetProperty("name").GetString()!,
                value.GetProperty("member").GetString()!,
                value.GetRawText()
            };
    }

    private static IDistribution Create(string name, JsonElement value)
    {
        var type = typeof(IDistribution).Assembly.GetType("UMapx.Distribution." + name)!;
        var signature = value.GetProperty("signature").EnumerateArray().Select(v => v.GetString() == "i" ? typeof(int) : typeof(float)).ToArray();
        var arguments = value.GetProperty("constructor").EnumerateArray().Select((v, i) => signature[i] == typeof(int) ? (object)(int)v.GetSingle() : v.GetSingle()).ToArray();
        return (IDistribution)type.GetConstructor(signature)!.Invoke(arguments);
    }

    private static double Evaluate(IDistribution instance, string member, JsonElement value)
    {
        var type = instance.GetType();
        try
        {
            return member switch
            {
                "Support.Min" => instance.Support.Min,
                "Support.Max" => instance.Support.Max,
                "Function" or "Distribution" => Convert.ToDouble(type.GetMethod(member, new[] { typeof(float) })!.Invoke(instance, new object[] { value.GetProperty("x").GetSingle() })),
                _ => Convert.ToDouble(type.GetProperty(member)!.GetValue(instance))
            };
        }
        catch (TargetInvocationException exception) when (exception.InnerException != null)
        {
            ExceptionDispatchInfo.Capture(exception.InnerException).Throw();
            throw;
        }
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void MatchesIndependentProbabilityReference(string distribution, string member, string json)
    {
        using var document = JsonDocument.Parse(json);
        var value = document.RootElement;
        var instance = Create(distribution, value);
        double actual = Evaluate(instance, member, value);
        var expectedElement = value.GetProperty("expected");
        if (expectedElement.ValueKind == JsonValueKind.String)
        {
            switch (expectedElement.GetString())
            {
                case "NaN":
                    Assert.True(double.IsNaN(actual), $"Expected undefined value; actual {actual}.");
                    return;
                case "Infinity":
                    Assert.True(double.IsPositiveInfinity(actual), $"Expected positive infinity; actual {actual}.");
                    return;
                case "-Infinity":
                    Assert.True(double.IsNegativeInfinity(actual), $"Expected negative infinity; actual {actual}.");
                    return;
            }
        }

        double expected = expectedElement.GetDouble(), tolerance = value.GetProperty("tolerance").GetDouble();
        // Preserve relative accuracy in positive tails, including rounding at the subnormal boundary.
        // Moment references can contain quadrature noise around an exact zero.
        double absolute = (member is "Function" or "Distribution") && expected > 0 && expected <= tolerance / (1 - tolerance) ? .5 * float.Epsilon : tolerance;
        Close(expected, actual, absolute, tolerance);
        if (instance is TukeyLambda tukey && member == "Function")
            Close(expected, tukey.Function(-value.GetProperty("x").GetSingle()), absolute, tolerance);
    }

    [Theory]
    [MemberData(nameof(EdgeCases))]
    public void MembersMatchHighPrecisionProbabilityReferences(string name, string member, string json)
    {
        using var document = JsonDocument.Parse(json);
        var value = document.RootElement;
        double actual = Evaluate(Create(name, value), member, value);
        var reference = value.GetProperty("expected");
        double expected = reference.ValueKind == JsonValueKind.String ? double.Parse(reference.GetString()!, System.Globalization.CultureInfo.InvariantCulture) : reference.GetDouble();
        if (double.IsNaN(expected))
            Assert.True(double.IsNaN(actual));
        else if (float.IsInfinity((float)expected))
            Assert.Equal((double)(float)expected, actual);
        else if (member == "Median" && name is "Poisson" or "Binomial")
            Assert.Equal(expected, actual);
        // Retain relative accuracy in small probabilities instead of accepting zero under a large absolute budget.
        else
            Close(expected, actual, 2 * (double)float.Epsilon, 2e-5);
    }
}
