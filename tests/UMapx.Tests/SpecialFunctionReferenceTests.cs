using System.Numerics;
using System.Text.Json;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Reference")]
public class SpecialFunctionReferenceTests
{
    public static IEnumerable<object[]> Cases()
    {
        foreach(var resource in new[]{"special-functions.json","special-functions-extended.json"})
        {
        using var stream = typeof(SpecialFunctionReferenceTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data."+resource);
        using var data = JsonDocument.Parse(stream!);
        foreach (var item in data.RootElement.EnumerateArray())
        {
            string name = item.GetProperty("name").GetString()!;
            string kinds = string.Join(",", item.GetProperty("kinds").EnumerateArray().Select(k => k.GetString()));
            yield return new object[] { name, kinds, item.GetProperty("args").GetRawText(), item.GetProperty("re").GetDouble(), item.GetProperty("im").GetDouble() };
        }
        }
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void AgreesWithHighPrecisionReference(string name, string signature, string arguments, double expectedReal, double expectedImaginary)
    {
        string[] kinds = signature.Split(',');
        using var data = JsonDocument.Parse(arguments);
        Type[] types = kinds.Select(k => k == "c" ? typeof(Complex32) : k == "i" ? typeof(int) : k == "b" ? typeof(bool) : typeof(float)).ToArray();
        object[] args = data.RootElement.EnumerateArray().Select((v, i) => kinds[i] == "c"
            ? (object)new Complex32(v[0].GetSingle(), v[1].GetSingle())
            : kinds[i] == "i" ? (object)v[0].GetInt32() : kinds[i] == "b" ? (object)(v[0].GetInt32()!=0) : v[0].GetSingle()).ToArray();
        var method = typeof(Special).GetMethod(name, types);
        Assert.NotNull(method);
        object? result = method.Invoke(null, args);
        Complex32 actual = result is Complex32 z ? z : new Complex32(Convert.ToSingle(result), 0);
        // A mixed absolute/relative bound, deliberately looser than float rounding.
        NumericAssert.Close(new Complex(expectedReal, expectedImaginary), actual, 2e-4, 2e-4);
    }
}
