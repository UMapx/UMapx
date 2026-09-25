using System.Numerics;
using System.Text.Json;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class SpecialFunctionReferenceTests
{
    public static IEnumerable<object[]> Cases() => ReadCases("special-functions.json", "special-functions-extended.json");
    public static IEnumerable<object[]> EdgeCases() => ReadCases("special-functions-edge-cases.json");
    private static IEnumerable<object[]> ReadCases(params string[] resources)
    {
        foreach (string resource in resources)
        {
            using var stream = typeof(SpecialFunctionReferenceTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data." + resource);
            using var data = JsonDocument.Parse(stream!);
            foreach (var item in data.RootElement.EnumerateArray())
                yield return new object[]
                {
                    item.GetProperty("name").GetString()!,
                    string.Join(",", item.GetProperty("kinds").EnumerateArray().Select(k => k.GetString())),
                    item.GetProperty("args").GetRawText(),
                    item.GetProperty("re").GetDouble(),
                    item.GetProperty("im").GetDouble()
                };
        }
    }

    private static object Evaluate(string name, string signature, string arguments)
    {
        string[] kinds = signature.Split(',');
        using var data = JsonDocument.Parse(arguments);
        Type[] types = kinds.Select(k => k == "c" ? typeof(Complex32) : k == "i" ? typeof(int) : k == "b" ? typeof(bool) : typeof(float)).ToArray();
        object[] args = data.RootElement.EnumerateArray().Select((v, i) => kinds[i] == "c" ? (object)new Complex32(v[0].GetSingle(), v[1].GetSingle()) : kinds[i] == "i" ? (object)v[0].GetInt32() : kinds[i] == "b" ? (object)(v[0].GetInt32() != 0) : v[0].GetSingle()).ToArray();
        var method = typeof(Special).GetMethod(name, types);
        Assert.NotNull(method);
        return method.Invoke(null, args)!;
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void AgreesWithHighPrecisionReference(string name, string signature, string arguments, double expectedReal, double expectedImaginary)
    {
        object result = Evaluate(name, signature, arguments);
        Complex32 actual = result is Complex32 z ? z : new Complex32(Convert.ToSingle(result), 0);
        var expected = new Complex(expectedReal, expectedImaginary);
        // Small representable values need relative accuracy; half a subnormal ULP permits rounding, not zeroing.
        // Keep the absolute budget at analytic zeros represented by negligible reference noise.
        double absolute = ((float)expectedReal != 0 || (float)expectedImaginary != 0) && expected.Magnitude <= 2e-4 / (1 - 2e-4) ? .5 * float.Epsilon : 2e-4;
        NumericAssert.Close(expected, actual, absolute, 2e-4);
    }

    [Theory]
    [MemberData(nameof(EdgeCases))]
    public void EdgeCasesMatchHighPrecisionReferences(string name, string signature, string arguments, double re, double im)
    {
        object result = Evaluate(name, signature, arguments);
        Complex actual = result is Complex32 z ? (Complex)z : new Complex(Convert.ToDouble(result), 0);
        var expected = new Complex(re, im);
        // Tiny tails must remain accurate; a fixed absolute tolerance would allow zero.
        double absolute = expected.Magnitude < 1e-4 && expected != Complex.Zero ? 1e-44 : 1e-7;
        double error = (actual - expected).Magnitude;
        Assert.True(double.IsFinite(error) && error <= absolute + 2e-5 * expected.Magnitude, $"Expected {expected}; actual {actual}; error {error:G9}.");
    }
}
