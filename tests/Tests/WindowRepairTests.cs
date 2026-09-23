using UMapx.Window;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Window")]
public class WindowRepairTests
{
    public static IEnumerable<object[]> Cases()
    {
        foreach (int size in new[] { 1, 2, 3, 8, 9, 32, 65 })
        foreach (string name in new[] { "Normal", "Confined", "BartlettHann" })
        foreach (float scale in new[] { .01f, .14f, 1f })
            yield return new object[] { name, size, scale };
    }

    [Theory, MemberData(nameof(Cases))]
    public void WindowsUseTheRequestedCenterAndRemainFiniteAtNarrowWidths(string name, int size, float scale)
    {
        WindowBase Create(int n) => name switch
        {
            "Normal" => new Normal(n, scale, 1.5f),
            "Confined" => new Confined(n, scale * size),
            _ => new BartlettHann(n)
        };
        var window = Create(size);
        var explicitSize = Create(17).GetWindow(size);
        var actual = window.GetWindow();
        Assert.Equal(size, actual.Length);
        Close(actual, explicitSize, 1e-7);
        for (int i = 0; i < size; i++)
        {
            Assert.True(float.IsFinite(actual[i]));
            Close(actual[i], actual[size - 1 - i], 2e-7, 0);
            double expected;
            if (size == 1) expected = 1;
            else if (name == "Normal")
                expected = Math.Exp(-Math.Pow(Math.Abs((i - (size - 1) / 2.0) / (scale * (size - 1) / 2.0)), 1.5));
            else if (name == "BartlettHann")
                expected = .62 - .48 * Math.Abs(i / (size - 1.0) - .5) - .38 * Math.Cos(2 * Math.PI * i / (size - 1));
            else
            {
                double sigma = scale * size;
                double G(double x) => Math.Exp(-Math.Pow((x - (size - 1) / 2.0) / (2 * sigma), 2));
                double edge = G(-.5);
                // At tiny widths both boundary samples underflow; their ratio tends to one.
                double ratio = edge == 0 ? 1 : edge / (G(size - .5) + G(-size - .5));
                expected = G(i) - ratio * (G(i + size) + G(i - size));
            }
            Close(expected, actual[i], 4e-7, 4e-7);
        }
    }
}
