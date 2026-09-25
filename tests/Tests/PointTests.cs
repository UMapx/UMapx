using System.Globalization;
using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class PointTests
{
    private static void Same((int X, int Y) expected, PointInt actual)
    {
        Assert.Equal(expected.X, actual.X);
        Assert.Equal(expected.Y, actual.Y);
    }

    private static void Same((float X, float Y) expected, PointFloat actual)
    {
        Assert.Equal(expected.X, actual.X);
        Assert.Equal(expected.Y, actual.Y);
    }

    private static readonly int[] Integers = { int.MinValue, -65536, -32769, -1, 0, 1, 32768, 65535, int.MaxValue };
    private static readonly float[] Floats = { float.MinValue, -3.5f, -float.Epsilon, -0.0f, 0, float.Epsilon,
        2.25f, float.MaxValue, float.NaN, float.PositiveInfinity, float.NegativeInfinity };

    [Fact]
    public void IntegerConstructionPropertiesEqualityAndOffsetsPreserveCoordinates()
    {
        foreach (int x in Integers)
        foreach (int y in Integers)
        {
            var expected = (X: x, Y: y);
            var actual = new PointInt(x, y);
            Same(expected, actual);
            Same(expected, new PointInt(new SizeInt(x, y)));
            Assert.Equal(x == 0 && y == 0, actual.IsEmpty);
            Assert.True(((IEquatable<PointInt>)actual).Equals(new PointInt(x, y)));
            Assert.True(actual.Equals((object)new PointInt(x, y)));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals((object)expected));
            Assert.Equal(x == 0 && y == 0, actual == PointInt.Empty);
            Assert.Equal(x != 0 || y != 0, actual != PointInt.Empty);
            Assert.Equal(new SizeInt(x, y), (SizeInt)actual);

            expected = (unchecked(x + y), unchecked(y + x));
            actual.Offset(new PointInt(y, x));
            Same(expected, actual);
            expected = (unchecked(expected.X - 10), unchecked(expected.Y + 17));
            actual.Offset(-10, 17);
            Same(expected, actual);
            actual.X = x; actual.Y = y;
            Same((x, y), actual);
        }
    }

    [Theory]
    [InlineData(0, 0, 0)]
    [InlineData(-1, -1, -1)]
    [InlineData(int.MinValue, 0, -32768)]
    [InlineData(int.MaxValue, -1, 32767)]
    [InlineData(65535, -1, 0)]
    [InlineData(65536, 0, 1)]
    [InlineData(-2147450880, -32768, -32768)]
    public void PackedIntegerConstructorDecodesSignedWords(int packed, int x, int y)
    {
        Same((x, y), new PointInt(packed));
    }

    [Fact]
    public void IntegerSizeArithmeticWrapsOverflowInEachCoordinate()
    {
        foreach (int x in Integers)
        foreach (int y in Integers)
        foreach (int delta in Integers)
        {
            var actual = new PointInt(x, y);
            var size = new SizeInt(delta, -7);
            var sum = (unchecked(x + delta), unchecked(y - 7));
            var difference = (unchecked(x - delta), unchecked(y + 7));
            Same(sum, PointInt.Add(actual, size));
            Same(difference, PointInt.Subtract(actual, size));
            Same(sum, actual + size);
            Same(difference, actual - size);
            Assert.Equal(new PointInt(x, y), actual);
        }
    }

    [Fact]
    public void FloatPropertiesEqualityAndSizeArithmeticUseEachCoordinate()
    {
        foreach (float x in Floats)
        foreach (float y in Floats)
        {
            var expected = (X: x, Y: y);
            var actual = new PointFloat(x, y);
            Same(expected, actual);
            Assert.Equal(x == 0 && y == 0, actual.IsEmpty);
            bool reflexive = !float.IsNaN(x) && !float.IsNaN(y);
            Assert.Equal(reflexive, actual.Equals(new PointFloat(x, y)));
            Assert.Equal(reflexive, actual.Equals((object)actual));
            Assert.Equal(reflexive, ((IEquatable<PointFloat>)actual).Equals(actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals((object)expected));
            Assert.Equal(x == 0 && y == 0, actual == PointFloat.Empty);
            Assert.Equal(x != 0 || y != 0, actual != PointFloat.Empty);

            var size = new SizeFloat(-2.75f, 4.5f);
            Same((x - 2.75f, y + 4.5f), PointFloat.Add(actual, size));
            Same((x + 2.75f, y - 4.5f), PointFloat.Subtract(actual, size));
            Same((x - 2.75f, y + 4.5f), actual + size);
            Same((x + 2.75f, y - 4.5f), actual - size);
            var integerSize = new SizeInt(int.MaxValue, int.MinValue);
            Same((x + int.MaxValue, y + int.MinValue), PointFloat.Add(actual, integerSize));
            Same((x - int.MaxValue, y - int.MinValue), PointFloat.Subtract(actual, integerSize));
            Same((x + int.MaxValue, y + int.MinValue), actual + integerSize);
            Same((x - int.MaxValue, y - int.MinValue), actual - integerSize);
            Same(expected, actual);
            actual.X = x; actual.Y = y;
            Same(expected, actual);
        }
        Assert.Equal(PointFloat.Empty.GetHashCode(), new PointFloat(-0.0f, -0.0f).GetHashCode());
    }

    [Theory]
    [InlineData(-2.5f, 3.5f, -2, 4, -2, 4, -2, 3)]
    [InlineData(2.5f, -3.5f, 3, -3, 2, -4, 2, -3)]
    [InlineData(-0.9f, 0.9f, 0, 1, -1, 1, 0, 0)]
    [InlineData(0, 16777215, 0, 16777215, 0, 16777215, 0, 16777215)]
    public void RoundingUsesCeilingMidpointToEvenAndTruncation(float x, float y,
        int ceilingX, int ceilingY, int roundX, int roundY, int truncateX, int truncateY)
    {
        var actual = new PointFloat(x, y);
        Same((ceilingX, ceilingY), PointInt.Ceiling(actual));
        Same((roundX, roundY), PointInt.Round(actual));
        Same((truncateX, truncateY), PointInt.Truncate(actual));
    }

    [Fact]
    public void NativeTypeAndVectorConversionsPreserveCoordinates()
    {
        foreach (int x in Integers)
        foreach (int y in Integers)
        {
            PointFloat point = new PointInt(x, y);
            Same(((float)x, (float)y), point);
        }
        foreach (float x in Floats)
        foreach (float y in Floats)
        {
            var vector = new Vector2(x, y);
            var point = new PointFloat(x, y);
            Same((x, y), new PointFloat(vector));
            Same((x, y), (PointFloat)vector);
            Assert.Equal(vector, point.ToVector2());
            Assert.Equal(vector, (Vector2)point);
        }
    }

    [Theory]
    [InlineData("en-US", "{X=-3.25, Y=7.5}")]
    [InlineData("ru-RU", "{X=-3,25, Y=7,5}")]
    public void StringFormattingUsesCurrentCultureAndCloneIsRetained(string culture, string expected)
    {
        var integer = new PointInt(-3, 7);
        var floating = new PointFloat(-3.25f, 7.5f);
        Assert.Equal(integer, integer.Clone());
        Assert.Equal(integer, (PointInt)((ICloneable)integer).Clone());
        Assert.Equal(floating, floating.Clone());
        Assert.Equal(floating, (PointFloat)((ICloneable)floating).Clone());
        var previous = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo(culture);
            Assert.Equal("{X=-3,Y=7}", integer.ToString());
            Assert.Equal(expected, floating.ToString());
        }
        finally { CultureInfo.CurrentCulture = previous; }
    }
}
