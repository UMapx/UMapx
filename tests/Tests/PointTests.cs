using System.Drawing;
using System.Globalization;
using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class PointTests
{
    private static void Same(Point expected, PointInt actual)
    {
        Assert.Equal(expected.X, actual.X);
        Assert.Equal(expected.Y, actual.Y);
    }

    private static void Same(PointF expected, PointFloat actual)
    {
        Assert.Equal(expected.X, actual.X);
        Assert.Equal(expected.Y, actual.Y);
    }

    private static readonly int[] Integers = { int.MinValue, -65536, -32769, -1, 0, 1, 32768, 65535, int.MaxValue };
    private static readonly float[] Floats = { float.MinValue, -3.5f, -float.Epsilon, -0.0f, 0, float.Epsilon,
        2.25f, float.MaxValue, float.NaN, float.PositiveInfinity, float.NegativeInfinity };

    [Fact]
    public void IntegerConstructionPropertiesEqualityAndOffsetsMatchSystemDrawing()
    {
        foreach (int x in Integers)
        foreach (int y in Integers)
        {
            var expected = new Point(x, y);
            var actual = new PointInt(x, y);
            Same(expected, actual);
            Same(new Point(new Size(x, y)), new PointInt(new SizeInt(x, y)));
            Assert.Equal(expected.IsEmpty, actual.IsEmpty);
            Assert.True(((IEquatable<PointInt>)actual).Equals(new PointInt(x, y)));
            Assert.True(actual.Equals((object)new PointInt(x, y)));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals((object)expected));
            Assert.Equal(expected == Point.Empty, actual == PointInt.Empty);
            Assert.Equal(expected != Point.Empty, actual != PointInt.Empty);
            Assert.Equal(new SizeInt(x, y), (SizeInt)actual);
            Assert.Equal(expected.GetHashCode(), actual.GetHashCode());
            Assert.Equal(expected.ToString(), actual.ToString());

            expected.Offset(new Point(y, x));
            actual.Offset(new PointInt(y, x));
            Same(expected, actual);
            expected.Offset(-10, 17);
            actual.Offset(-10, 17);
            Same(expected, actual);
            actual.X = x; actual.Y = y;
            Same(new Point(x, y), actual);
        }
    }

    [Theory]
    [InlineData(0)]
    [InlineData(-1)]
    [InlineData(int.MinValue)]
    [InlineData(int.MaxValue)]
    [InlineData(65535)]
    [InlineData(65536)]
    [InlineData(-2147450880)]
    public void PackedIntegerConstructorMatchesSignedWordsInSystemDrawing(int packed)
    {
        Same(new Point(packed), new PointInt(packed));
    }

    [Fact]
    public void IntegerSizeArithmeticMatchesSystemDrawingIncludingOverflow()
    {
        foreach (int x in Integers)
        foreach (int y in Integers)
        foreach (int delta in Integers)
        {
            var expected = new Point(x, y);
            var actual = new PointInt(x, y);
            var drawingSize = new Size(delta, -7);
            var size = new SizeInt(delta, -7);
            Same(Point.Add(expected, drawingSize), PointInt.Add(actual, size));
            Same(Point.Subtract(expected, drawingSize), PointInt.Subtract(actual, size));
            Same(expected + drawingSize, actual + size);
            Same(expected - drawingSize, actual - size);
            Assert.Equal(new PointInt(x, y), actual);
        }
    }

    [Fact]
    public void FloatPropertiesEqualityAndSizeArithmeticMatchSystemDrawing()
    {
        foreach (float x in Floats)
        foreach (float y in Floats)
        {
            var expected = new PointF(x, y);
            var actual = new PointFloat(x, y);
            Same(expected, actual);
            Assert.Equal(expected.IsEmpty, actual.IsEmpty);
            Assert.Equal(expected.Equals(new PointF(x, y)), actual.Equals(new PointFloat(x, y)));
            Assert.Equal(expected.Equals((object)expected), actual.Equals((object)actual));
            Assert.Equal(expected.Equals(expected), ((IEquatable<PointFloat>)actual).Equals(actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals((object)expected));
            Assert.Equal(expected == PointF.Empty, actual == PointFloat.Empty);
            Assert.Equal(expected != PointF.Empty, actual != PointFloat.Empty);
            Assert.Equal(expected.GetHashCode(), actual.GetHashCode());
            Assert.Equal(expected.ToString(), actual.ToString());

            var size = new SizeFloat(-2.75f, 4.5f);
            var drawingSize = new SizeF(size.Width, size.Height);
            Same(PointF.Add(expected, drawingSize), PointFloat.Add(actual, size));
            Same(PointF.Subtract(expected, drawingSize), PointFloat.Subtract(actual, size));
            Same(expected + drawingSize, actual + size);
            Same(expected - drawingSize, actual - size);
            var integerSize = new SizeInt(int.MaxValue, int.MinValue);
            var drawingIntegerSize = new Size(integerSize.Width, integerSize.Height);
            Same(PointF.Add(expected, drawingIntegerSize), PointFloat.Add(actual, integerSize));
            Same(PointF.Subtract(expected, drawingIntegerSize), PointFloat.Subtract(actual, integerSize));
            Same(expected + drawingIntegerSize, actual + integerSize);
            Same(expected - drawingIntegerSize, actual - integerSize);
            Same(expected, actual);
            actual.X = x; actual.Y = y;
            Same(expected, actual);
        }
        Assert.Equal(PointFloat.Empty.GetHashCode(), new PointFloat(-0.0f, -0.0f).GetHashCode());
    }

    [Theory]
    [InlineData(-2.5f, 3.5f)]
    [InlineData(2.5f, -3.5f)]
    [InlineData(-0.9f, 0.9f)]
    [InlineData(0, 16777215)]
    [InlineData(float.NaN, float.PositiveInfinity)]
    [InlineData(float.MaxValue, float.MinValue)]
    public void RoundingMatchesSystemDrawing(float x, float y)
    {
        var expected = new PointF(x, y);
        var actual = new PointFloat(x, y);
        Same(Point.Ceiling(expected), PointInt.Ceiling(actual));
        Same(Point.Round(expected), PointInt.Round(actual));
        Same(Point.Truncate(expected), PointInt.Truncate(actual));
    }

    [Fact]
    public void NativeTypeAndVectorConversionsPreserveCoordinates()
    {
        foreach (int x in Integers)
        foreach (int y in Integers)
        {
            PointFloat point = new PointInt(x, y);
            Same((PointF)new Point(x, y), point);
        }
        foreach (float x in Floats)
        foreach (float y in Floats)
        {
            var vector = new Vector2(x, y);
            var point = new PointFloat(x, y);
            Same(new PointF(x, y), new PointFloat(vector));
            Same(new PointF(x, y), (PointFloat)vector);
            Assert.Equal(vector, point.ToVector2());
            Assert.Equal(vector, (Vector2)point);
        }
    }

    [Theory]
    [InlineData("en-US")]
    [InlineData("ru-RU")]
    public void StringFormattingMatchesSystemDrawingAndCloneIsRetained(string culture)
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
            Assert.Equal(new Point(-3, 7).ToString(), integer.ToString());
            Assert.Equal(new PointF(-3.25f, 7.5f).ToString(), floating.ToString());
        }
        finally { CultureInfo.CurrentCulture = previous; }
    }
}
