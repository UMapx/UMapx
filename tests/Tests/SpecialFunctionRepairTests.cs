using System.Numerics;
using System.Text.Json;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Reference")]
public class SpecialFunctionRepairTests
{
    public static IEnumerable<object[]> Cases()
    {
        using var stream = typeof(SpecialFunctionRepairTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.special-functions-repair.json");
        using var data = JsonDocument.Parse(stream!);
        foreach (var item in data.RootElement.EnumerateArray())
            yield return new object[] { item.GetProperty("name").GetString()!,
                string.Join(",", item.GetProperty("kinds").EnumerateArray().Select(k => k.GetString())),
                item.GetProperty("args").GetRawText(), item.GetProperty("re").GetDouble(), item.GetProperty("im").GetDouble() };
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void AdditionalHighPrecisionReferences(string name, string signature, string arguments, double re, double im)
    {
        var kinds = signature.Split(',');
        using var data = JsonDocument.Parse(arguments);
        var types = kinds.Select(k => k == "c" ? typeof(Complex32) : k == "i" ? typeof(int) : k == "b" ? typeof(bool) : typeof(float)).ToArray();
        var args = data.RootElement.EnumerateArray().Select((v, i) => kinds[i] == "c" ? (object)new Complex32(v[0].GetSingle(), v[1].GetSingle())
            : kinds[i] == "i" ? (object)v[0].GetInt32() : kinds[i] == "b" ? (object)(v[0].GetInt32() != 0) : v[0].GetSingle()).ToArray();
        var method = typeof(Special).GetMethod(name, types);
        Assert.NotNull(method);
        var value = method.Invoke(null, args);
        Complex actual = value is Complex32 z ? (Complex)z : new Complex(Convert.ToDouble(value), 0);
        var expected = new Complex(re, im);
        // Tiny tails must remain accurate; a fixed absolute tolerance would allow zero.
        double absolute = expected.Magnitude < 1e-4 && expected != Complex.Zero ? 1e-44 : 1e-7;
        double error = (actual - expected).Magnitude;
        Assert.True(double.IsFinite(error) && error <= absolute + 2e-5 * expected.Magnitude,
            $"Expected {expected}; actual {actual}; error {error:G9}.");
    }
}
