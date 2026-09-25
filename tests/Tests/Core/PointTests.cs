using System.Globalization;
using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Geometry")]
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

    private static readonly int[] Integers =
    {
        int.MinValue,
        -65536,
        -32769,
        -1,
        0,
        1,
        32768,
        65535,
        int.MaxValue
    };
    private static readonly float[] Floats =
    {
        float.MinValue,
        -3.5f,
        -float.Epsilon,
        -0.0f,
        0,
        float.Epsilon,
        2.25f,
        float.MaxValue,
        float.NaN,
        float.PositiveInfinity,
        float.NegativeInfinity
    };
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
                actual.X = x;
                actual.Y = y;
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
                actual.X = x;
                actual.Y = y;
                Same(expected, actual);
            }

        Assert.Equal(PointFloat.Empty.GetHashCode(), new PointFloat(-0.0f, -0.0f).GetHashCode());
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
        finally
        {
            CultureInfo.CurrentCulture = previous;
        }
    }

    [Fact]
    public void IntegerTranslationReturnsNewArraysAndPreservesOrderAndInput()
    {
        var points = new[]
        {
            new PointInt(-3, 5),
            new PointInt(2, -2)
        };
        var original = (PointInt[])points.Clone();
        var offset = new PointInt(4, -1);
        var added = PointInt.Add(points, offset);
        var subtracted = PointInt.Sub(points, offset);
        Assert.Equal(new[] { new PointInt(1, 4), new PointInt(6, -3) }, added);
        Assert.Equal(new[] { new PointInt(-7, 6), new PointInt(-2, -1) }, subtracted);
        Assert.NotSame(points, added);
        Assert.NotSame(points, subtracted);
        Assert.Equal(original, points);
    }

    [Fact]
    public void FloatTranslationPreservesFractionalCoordinatesAndInput()
    {
        var points = new[]
        {
            new PointFloat(-3.25f, 5.5f),
            new PointFloat(2.5f, -2.25f)
        };
        var original = (PointFloat[])points.Clone();
        var offset = new PointFloat(4.5f, -1.25f);
        var added = PointFloat.Add(points, offset);
        var subtracted = PointFloat.Sub(points, offset);
        Assert.Equal(new[] { new PointFloat(1.25f, 4.25f), new PointFloat(7, -3.5f) }, added);
        Assert.Equal(new[] { new PointFloat(-7.75f, 6.75f), new PointFloat(-2, -1) }, subtracted);
        Assert.NotSame(points, added);
        Assert.NotSame(points, subtracted);
        Assert.Equal(original, points);
    }

    [Fact]
    public void RotationUsesDegreesAndTheSpecifiedCenter()
    {
        var center = new PointInt(1, 1);
        var points = new[]
        {
            new PointInt(3, 2),
            center
        };
        var rotated = PointInt.Rotate(points, center, 90);
        Assert.Equal(new[] { new PointInt(0, 3), center }, rotated);
        Assert.Equal(new PointInt(0, 3), points[0].Rotate(center, 90));
        Assert.Equal(new PointInt(1, 1), new PointInt(2, 0).Rotate(PointInt.Empty, 60));
        Assert.NotSame(points, rotated);
        Assert.Equal(new[] { new PointInt(3, 2), center }, points);
        var floatingCenter = new PointFloat(1.25f, 1.5f);
        var floatingPoints = new[]
        {
            new PointFloat(3.25f, 2.5f),
            floatingCenter
        };
        var floatingRotated = PointFloat.Rotate(floatingPoints, floatingCenter, 90);
        Assert.Equal(new[] { new PointFloat(0.25f, 3.5f), floatingCenter }, floatingRotated);
        Assert.Equal(new PointFloat(0.25f, 3.5f), floatingPoints[0].Rotate(floatingCenter, 90));
        var sixtyDegrees = new PointFloat(2, 0).Rotate(PointFloat.Empty, 60);
        NumericAssert.Close(1, sixtyDegrees.X, 1e-6, 0);
        NumericAssert.Close(1.7320508075688772, sixtyDegrees.Y, 1e-6, 0);
        Assert.NotSame(floatingPoints, floatingRotated);
        Assert.Equal(new[] { new PointFloat(3.25f, 2.5f), floatingCenter }, floatingPoints);
    }

    [Fact]
    public void BoundsUseCoordinateExtremaAndMeanTruncatesOnlyForIntegers()
    {
        var points = new[]
        {
            new PointInt(-3, 5),
            new PointInt(2, -2)
        };
        Assert.Equal(new RectangleInt(-3, -2, 5, 7), PointInt.GetRectangle(points));
        Assert.Equal(new PointInt(0, 1), PointInt.GetMeanPoint(points));
        Assert.Equal(new PointInt(2, 5), points[0].GetSupportedPoint(points[1]));
        var floating = new[]
        {
            new PointFloat(-3.25f, 5.5f),
            new PointFloat(2.5f, -2.25f)
        };
        Assert.Equal(new RectangleFloat(-3.25f, -2.25f, 5.75f, 7.75f), PointFloat.GetRectangle(floating));
        Assert.Equal(new PointFloat(-0.375f, 1.625f), PointFloat.GetMeanPoint(floating));
        Assert.Equal(new PointFloat(2.5f, 5.5f), floating[0].GetSupportedPoint(floating[1]));
    }

    [Theory]
    [InlineData(1, -89.99337047)]
    [InlineData(-1, 89.99337047)]
    public void AnglePreservesSignAndApproximateDegreeConversion(int rightY, double expected)
    {
        NumericAssert.Close(expected, PointInt.Empty.GetAngle(new PointInt(0, rightY), new PointInt(1, 0)), 1e-5, 0);
        NumericAssert.Close(expected, PointFloat.Empty.GetAngle(new PointFloat(0, rightY), new PointFloat(1, 0)), 1e-5, 0);
    }

    [Fact]
    public void DegenerateAnglesPreserveDivisionBehavior()
    {
        Assert.True(float.IsNaN(PointInt.Empty.GetAngle(PointInt.Empty, PointInt.Empty)));
        Assert.True(float.IsNaN(PointFloat.Empty.GetAngle(PointFloat.Empty, PointFloat.Empty)));
        NumericAssert.Close(-89.99337047, PointInt.Empty.GetAngle(new PointInt(1, 0), PointInt.Empty), 1e-5, 0);
        NumericAssert.Close(-89.99337047, PointFloat.Empty.GetAngle(new PointFloat(1, 0), PointFloat.Empty), 1e-5, 0);
    }

    [Fact]
    public void EmptyArraysPreserveArithmeticResults()
    {
        Assert.Empty(PointInt.Add(Array.Empty<PointInt>(), new PointInt(1, 2)));
        Assert.Empty(PointInt.Sub(Array.Empty<PointInt>(), new PointInt(1, 2)));
        Assert.Empty(PointInt.Rotate(Array.Empty<PointInt>(), PointInt.Empty, 90));
        Assert.Equal(new RectangleInt(int.MaxValue, int.MaxValue, 1, 1), PointInt.GetRectangle(Array.Empty<PointInt>()));
        Assert.Throws<DivideByZeroException>(() => PointInt.GetMeanPoint());
        Assert.Empty(PointFloat.Add(Array.Empty<PointFloat>(), new PointFloat(1, 2)));
        Assert.Empty(PointFloat.Sub(Array.Empty<PointFloat>(), new PointFloat(1, 2)));
        Assert.Empty(PointFloat.Rotate(Array.Empty<PointFloat>(), PointFloat.Empty, 90));
        Assert.Equal(new RectangleFloat(float.MaxValue, float.MaxValue, float.NegativeInfinity, float.NegativeInfinity), PointFloat.GetRectangle(Array.Empty<PointFloat>()));
        var mean = PointFloat.GetMeanPoint();
        Assert.True(float.IsNaN(mean.X));
        Assert.True(float.IsNaN(mean.Y));
    }

    [Fact]
    public void IntegerArithmeticPreservesOverflow()
    {
        var point = new PointInt(int.MaxValue, int.MinValue);
        Assert.Equal(new[] { new PointInt(int.MinValue, int.MaxValue) }, PointInt.Add(new[] { point }, new PointInt(1, -1)));
        Assert.Equal(new[] { new PointInt(int.MinValue, int.MaxValue) }, PointInt.Sub(new[] { point }, new PointInt(-1, 1)));
        Assert.Equal(new PointInt(-1, 0), PointInt.GetMeanPoint(point, point));
        Assert.Equal(new RectangleInt(int.MinValue, int.MinValue, -1, -1), PointInt.GetRectangle(new[] { point, new PointInt(int.MinValue, int.MaxValue) }));
    }
}
