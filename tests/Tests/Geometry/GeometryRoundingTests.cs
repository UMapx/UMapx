using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Geometry")]
public class GeometryRoundingTests
{
    [Theory]
    [InlineData(-3.5f, -3, -4, -3)]
    [InlineData(-2.5f, -2, -2, -2)]
    [InlineData(-0.9f, 0, -1, 0)]
    [InlineData(-float.Epsilon, 0, 0, 0)]
    [InlineData(0f, 0, 0, 0)]
    [InlineData(float.Epsilon, 1, 0, 0)]
    [InlineData(0.9f, 1, 1, 0)]
    [InlineData(2.5f, 3, 2, 2)]
    [InlineData(3.5f, 4, 4, 3)]
    [InlineData(16777215f, 16777215, 16777215, 16777215)]
    [InlineData(3.75f, 4, 4, 3)]
    [InlineData(-1.25f, -1, -1, -1)]
    [InlineData(-2.9f, -2, -3, -2)]
    [InlineData(3.1f, 4, 3, 3)]
    [InlineData(0.1f, 1, 0, 0)]
    [InlineData(-0.1f, 0, 0, 0)]
    [InlineData(4.5f, 5, 4, 4)]
    [InlineData(-5.5f, -5, -6, -5)]
    [InlineData(100000f, 100000, 100000, 100000)]
    [InlineData(200000f, 200000, 200000, 200000)]
    [InlineData(-16777215f, -16777215, -16777215, -16777215)]
    public void RoundingUsesSpecifiedRulesInEveryComponent(float value, int ceiling, int round, int truncate)
    {
        Check(value, ceiling, round, truncate);
    }

    [Theory]
    [InlineData(float.NaN)]
    [InlineData(float.PositiveInfinity)]
    [InlineData(float.NegativeInfinity)]
    [InlineData(float.MinValue)]
    [InlineData(float.MaxValue)]
    [InlineData(2147483648f)]
    [InlineData(-2147483648f)]
    public void ExceptionalComponentsUseUncheckedRuntimeNumericConversions(float value)
    {
        Check(value, unchecked((int)Math.Ceiling(value)), unchecked((int)Math.Round(value)), unchecked((int)value));
    }

    private static void Check(float value, int ceiling, int round, int truncate)
    {
        Assert.Equal(new PointInt(ceiling, 1), PointInt.Ceiling(new PointFloat(value, 1)));
        Assert.Equal(new PointInt(-1, ceiling), PointInt.Ceiling(new PointFloat(-1, value)));
        Assert.Equal(new PointInt(round, 1), PointInt.Round(new PointFloat(value, 1)));
        Assert.Equal(new PointInt(-1, round), PointInt.Round(new PointFloat(-1, value)));
        Assert.Equal(new PointInt(truncate, 1), PointInt.Truncate(new PointFloat(value, 1)));
        Assert.Equal(new PointInt(-1, truncate), PointInt.Truncate(new PointFloat(-1, value)));
        Assert.Equal(new SizeInt(ceiling, 1), SizeInt.Ceiling(new SizeFloat(value, 1)));
        Assert.Equal(new SizeInt(-1, ceiling), SizeInt.Ceiling(new SizeFloat(-1, value)));
        Assert.Equal(new SizeInt(round, 1), SizeInt.Round(new SizeFloat(value, 1)));
        Assert.Equal(new SizeInt(-1, round), SizeInt.Round(new SizeFloat(-1, value)));
        Assert.Equal(new SizeInt(truncate, 1), SizeInt.Truncate(new SizeFloat(value, 1)));
        Assert.Equal(new SizeInt(-1, truncate), SizeInt.Truncate(new SizeFloat(-1, value)));
        Assert.Equal(new SizeInt(truncate, 1), new SizeFloat(value, 1).ToSize());
        Assert.Equal(new SizeInt(-1, truncate), new SizeFloat(-1, value).ToSize());
        Assert.Equal(new RangeInt(ceiling, 1), RangeInt.Ceiling(new RangeFloat(value, 1)));
        Assert.Equal(new RangeInt(-1, ceiling), RangeInt.Ceiling(new RangeFloat(-1, value)));
        Assert.Equal(new RangeInt(round, 1), RangeInt.Round(new RangeFloat(value, 1)));
        Assert.Equal(new RangeInt(-1, round), RangeInt.Round(new RangeFloat(-1, value)));
        Assert.Equal(new RangeInt(truncate, 1), RangeInt.Truncate(new RangeFloat(value, 1)));
        Assert.Equal(new RangeInt(-1, truncate), RangeInt.Truncate(new RangeFloat(-1, value)));
        for (int component = 0; component < 4; component++)
        {
            float[] values =
            {
                -1,
                2,
                3,
                4
            };
            values[component] = value;
            var rectangle = new RectangleFloat(values[0], values[1], values[2], values[3]);
            foreach (var (expected, actual) in new[]
            {
                (ceiling, RectangleInt.Ceiling(rectangle)),
                (round, RectangleInt.Round(rectangle)),
                (truncate, RectangleInt.Truncate(rectangle))
            }

            )
            {
                int[] components =
                {
                    -1,
                    2,
                    3,
                    4
                };
                components[component] = expected;
                Assert.Equal(new RectangleInt(components[0], components[1], components[2], components[3]), actual);
            }
        }
    }
}
