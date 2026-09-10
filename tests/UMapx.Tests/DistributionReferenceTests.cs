using System.Reflection;
using System.Text.Json;
using UMapx.Core;
using UMapx.Distribution;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Distribution")]
public class DistributionReferenceTests
{
    public static IEnumerable<object[]> Cases()
    {
        using var stream = typeof(DistributionReferenceTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.distributions.json");
        using var document = JsonDocument.Parse(stream!);
        foreach (var value in document.RootElement.EnumerateArray())
            yield return new object[] { value.GetProperty("name").GetString()!, value.GetProperty("member").GetString()!, value.GetRawText() };
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void MatchesIndependentProbabilityReference(string distribution, string member, string json)
    {
        using var document = JsonDocument.Parse(json);
        var value = document.RootElement;
        var type = typeof(IDistribution).Assembly.GetType("UMapx.Distribution." + distribution)!;
        var signature = value.GetProperty("signature").EnumerateArray().Select(v => v.GetString() == "i" ? typeof(int) : typeof(float)).ToArray();
        var arguments = value.GetProperty("constructor").EnumerateArray().Select((v, i) => signature[i] == typeof(int) ? (object)(int)v.GetSingle() : v.GetSingle()).ToArray();
        var instance = (IDistribution)type.GetConstructor(signature)!.Invoke(arguments);
        double actual;
        try
        {
            actual = member switch
            {
                "Support.Min" => instance.Support.Min,
                "Support.Max" => instance.Support.Max,
                "Function" or "Distribution" => Convert.ToDouble(type.GetMethod(member, new[] { typeof(float) })!.Invoke(instance, new object[] { value.GetProperty("x").GetSingle() })),
                _ => Convert.ToDouble(type.GetProperty(member)!.GetValue(instance))
            };
        }
        catch (TargetInvocationException exception) when (exception.InnerException != null)
        {
            throw exception.InnerException;
        }
        var expectedElement = value.GetProperty("expected");
        if (expectedElement.ValueKind == JsonValueKind.String)
        {
            switch (expectedElement.GetString())
            {
                case "NaN": Assert.True(double.IsNaN(actual), $"Expected undefined value; actual {actual}."); return;
                case "Infinity": Assert.True(double.IsPositiveInfinity(actual), $"Expected positive infinity; actual {actual}."); return;
                case "-Infinity": Assert.True(double.IsNegativeInfinity(actual), $"Expected negative infinity; actual {actual}."); return;
            }
        }
        NumericAssert.Close(expectedElement.GetDouble(), actual, value.GetProperty("tolerance").GetDouble(), value.GetProperty("tolerance").GetDouble());
    }
}
