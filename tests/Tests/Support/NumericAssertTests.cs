using Xunit.Sdk;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Infrastructure")]
public class NumericAssertTests
{
    [Theory]
    [InlineData(double.PositiveInfinity, 0)]
    [InlineData(double.NegativeInfinity, 0)]
    [InlineData(double.NaN, 0)]
    [InlineData(0, double.PositiveInfinity)]
    [InlineData(0, double.NegativeInfinity)]
    [InlineData(0, double.NaN)]
    public void CloseRejectsNonFiniteValues(double expected, double actual)
    {
        Assert.ThrowsAny<XunitException>(() => NumericAssert.Close(expected, actual));
    }
}
