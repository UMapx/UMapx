using System.Drawing;
using System.Globalization;
using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class RectangleTests
{
    private static void Same(Rectangle expected, RectangleInt actual)
    {
        Assert.Equal(expected.X, actual.X);
        Assert.Equal(expected.Y, actual.Y);
        Assert.Equal(expected.Width, actual.Width);
        Assert.Equal(expected.Height, actual.Height);
    }

    private static void Same(RectangleF expected, RectangleFloat actual)
    {
        Assert.Equal(expected.X, actual.X);
        Assert.Equal(expected.Y, actual.Y);
        Assert.Equal(expected.Width, actual.Width);
        Assert.Equal(expected.Height, actual.Height);
    }

    private static Rectangle[] IntegerCases()
    {
        var cases = new List<Rectangle>
        {
            Rectangle.Empty, new(3, 4, 0, 0), new(0, 0, 10, 20), new(10, 20, 5, 5),
            new(10, 0, 4, 3), new(-4, -3, 10, 20), new(2, 3, -4, 5), new(2, 3, 4, -5),
            new(2, 3, -4, -5), new(5, 5, 0, 4), new(5, 5, 4, 0),
            new(int.MaxValue, int.MinValue, 4, -5), new(int.MinValue, int.MaxValue, -4, 5),
            new(-1, -2, int.MaxValue, int.MinValue)
        };
        var random = new Random(731);
        for (int i = 0; i < 60; i++)
            cases.Add(new Rectangle(random.Next(-50, 50), random.Next(-50, 50),
                random.Next(-10, 70), random.Next(-10, 70)));
        return cases.ToArray();
    }

    private static RectangleF[] FloatCases()
    {
        var cases = IntegerCases().Select(r => (RectangleF)r).ToList();
        cases.AddRange(new[]
        {
            new RectangleF(-2.75f, 3.125f, 4.5f, 1.75f),
            new RectangleF(-0.0f, 0, 0, -0.0f),
            new RectangleF(float.Epsilon, -float.Epsilon, float.Epsilon, float.Epsilon),
            new RectangleF(float.MaxValue, float.MinValue, float.MaxValue, float.MaxValue),
            new RectangleF(float.PositiveInfinity, 0, 2, 3),
            new RectangleF(0, float.NegativeInfinity, 2, 3),
            new RectangleF(0, 0, float.PositiveInfinity, 3),
            new RectangleF(float.NaN, 0, 2, 3),
            new RectangleF(0, float.NaN, 2, 3),
            new RectangleF(0, 0, float.NaN, 3),
            new RectangleF(0, 0, 2, float.NaN)
        });
        return cases.ToArray();
    }

    [Fact]
    public void IntegerConstructionPropertiesAndEqualityMatchSystemDrawing()
    {
        foreach (var expected in IntegerCases())
        {
            var actual = new RectangleInt(new PointInt(expected.X, expected.Y),
                new SizeInt(expected.Width, expected.Height));
            Same(expected, actual);
            Assert.Equal(expected.Left, actual.Left);
            Assert.Equal(expected.Top, actual.Top);
            Assert.Equal(expected.Right, actual.Right);
            Assert.Equal(expected.Bottom, actual.Bottom);
            Assert.Equal(expected.ToString(), actual.ToString());
            Assert.Equal(new PointInt(expected.X, expected.Y), actual.Location);
            Assert.Equal(new SizeInt(expected.Width, expected.Height), actual.Size);
            Same(Rectangle.FromLTRB(expected.X, expected.Y, expected.Width, expected.Height),
                RectangleInt.FromLTRB(expected.X, expected.Y, expected.Width, expected.Height));
            Assert.True(actual.Equals((object)actual));
            Assert.True(((IEquatable<RectangleInt>)actual).Equals(actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals((object)expected));
            Assert.Equal(actual.GetHashCode(), actual.Clone().GetHashCode());
            Same(expected, actual.Clone());
            Same(expected, (RectangleInt)((ICloneable)actual).Clone());

            var changed = actual;
            changed.Location = new PointInt(7, -11);
            changed.Size = new SizeInt(13, 17);
            Same(new Rectangle(7, -11, 13, 17), changed);
            changed.X = -9; changed.Y = 3; changed.Width = 0; changed.Height = -4;
            Same(new Rectangle(-9, 3, 0, -4), changed);
            Same(expected, actual);
        }
    }

    [Fact]
    public void FloatConstructionPropertiesAndEqualityMatchSystemDrawing()
    {
        foreach (var expected in FloatCases())
        {
            var actual = new RectangleFloat(new PointFloat(expected.X, expected.Y),
                new SizeFloat(expected.Width, expected.Height));
            Same(expected, actual);
            Assert.Equal(expected.Left, actual.Left);
            Assert.Equal(expected.Top, actual.Top);
            Assert.Equal(expected.Right, actual.Right);
            Assert.Equal(expected.Bottom, actual.Bottom);
            Assert.Equal(expected.IsEmpty, actual.IsEmpty);
            Assert.Equal(expected.ToString(), actual.ToString());
            Assert.Equal(expected.X, actual.Location.X);
            Assert.Equal(expected.Y, actual.Location.Y);
            Assert.Equal(expected.Width, actual.Size.Width);
            Assert.Equal(expected.Height, actual.Size.Height);
            Same(RectangleF.FromLTRB(expected.X, expected.Y, expected.Width, expected.Height),
                RectangleFloat.FromLTRB(expected.X, expected.Y, expected.Width, expected.Height));
            Assert.Equal(expected.Equals(expected), actual.Equals((object)actual));
            Assert.Equal(expected.Equals(expected), ((IEquatable<RectangleFloat>)actual).Equals(actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals((object)expected));
            Assert.Equal(actual.GetHashCode(), actual.Clone().GetHashCode());
            Same(expected, actual.Clone());
            Same(expected, (RectangleFloat)((ICloneable)actual).Clone());

            var changed = actual;
            changed.Location = new PointFloat(7.25f, -11.5f);
            changed.Size = new SizeFloat(13.75f, 17.5f);
            Same(new RectangleF(7.25f, -11.5f, 13.75f, 17.5f), changed);
            changed.X = -9.5f; changed.Y = 3.25f; changed.Width = 0; changed.Height = -4;
            Same(new RectangleF(-9.5f, 3.25f, 0, -4), changed);
            Same(expected, actual);
        }
        Assert.Equal(RectangleFloat.Empty.GetHashCode(), new RectangleFloat(-0.0f, 0, 0, -0.0f).GetHashCode());
    }

    [Theory]
    [InlineData(0, 0, 0, 0, true)]
    [InlineData(10, 20, 0, 0, true)]
    [InlineData(10, 20, 0, 5, true)]
    [InlineData(10, 20, 5, 0, true)]
    [InlineData(10, 20, -5, 7, true)]
    [InlineData(10, 20, 5, -7, true)]
    [InlineData(10, 20, -5, -7, true)]
    [InlineData(0, 0, 1, 1, false)]
    [InlineData(-10, -20, 5, 7, false)]
    [InlineData(int.MaxValue, int.MinValue, 5, 7, false)]
    public void BothTypesUseDimensionsToDetermineEmptiness(int x, int y, int width, int height, bool expected)
    {
        Assert.Equal(expected, new RectangleInt(x, y, width, height).IsEmpty);
        Assert.Equal(expected, new RectangleFloat(x, y, width, height).IsEmpty);
    }

    [Fact]
    public void TouchingIntersectionIsEmptyWithoutBeingTheDefaultValue()
    {
        var integer = RectangleInt.Intersect(new RectangleInt(0, 0, 10, 10), new RectangleInt(10, 0, 5, 10));
        Assert.True(integer.IsEmpty);
        Assert.NotEqual(RectangleInt.Empty, integer);
        Assert.Equal(new RectangleInt(10, 0, 0, 10), integer);

        var floating = RectangleFloat.Intersect(new RectangleFloat(0, 0, 10, 10), new RectangleFloat(10, 0, 5, 10));
        Assert.True(floating.IsEmpty);
        Assert.NotEqual(RectangleFloat.Empty, floating);
        Assert.Equal(new RectangleFloat(10, 0, 0, 10), floating);
    }

    [Fact]
    public void IntegerGeometryMatchesSystemDrawingIncludingEmptyNegativeAndOverflowingBounds()
    {
        var cases = IntegerCases();
        foreach (var a in cases)
        foreach (var b in cases)
        {
            RectangleInt actual = a;
            RectangleInt other = b;
            Assert.Equal(a == b, actual == other);
            Assert.Equal(a != b, actual != other);
            Assert.Equal(a.Equals(b), actual.Equals(other));
            Assert.Equal(a.Contains(b), actual.Contains(other));
            Assert.Equal(a.Contains(b.X, b.Y), actual.Contains(b.X, b.Y));
            Assert.Equal(a.Contains(b.Location), actual.Contains(new PointInt(b.X, b.Y)));
            Assert.Equal(a.IntersectsWith(b), actual.IntersectsWith(other));
            Same(Rectangle.Intersect(a, b), RectangleInt.Intersect(actual, other));
            Same(Rectangle.Union(a, b), RectangleInt.Union(actual, other));
            var expected = a;
            expected.Intersect(b);
            actual.Intersect(other);
            Same(expected, actual);
        }
    }

    [Fact]
    public void FloatGeometryMatchesSystemDrawingIncludingNaNAndInfinity()
    {
        var cases = FloatCases();
        foreach (var a in cases)
        foreach (var b in cases)
        {
            RectangleFloat actual = a;
            RectangleFloat other = b;
            Assert.Equal(a == b, actual == other);
            Assert.Equal(a != b, actual != other);
            Assert.Equal(a.Equals(b), actual.Equals(other));
            Assert.Equal(a.Contains(b), actual.Contains(other));
            Assert.Equal(a.Contains(b.X, b.Y), actual.Contains(b.X, b.Y));
            Assert.Equal(a.Contains(b.Location), actual.Contains(new PointFloat(b.X, b.Y)));
            Assert.Equal(a.IntersectsWith(b), actual.IntersectsWith(other));
            Same(RectangleF.Intersect(a, b), RectangleFloat.Intersect(actual, other));
            Same(RectangleF.Union(a, b), RectangleFloat.Union(actual, other));
            var expected = a;
            expected.Intersect(b);
            actual.Intersect(other);
            Same(expected, actual);
        }
    }

    [Fact]
    public void IntegerOffsetAndInflationMatchSystemDrawing()
    {
        foreach (var rectangle in IntegerCases())
        foreach (int amount in new[] { -20, 0, 7, int.MinValue, int.MaxValue })
        {
            RectangleInt original = rectangle;
            Same(Rectangle.Inflate(rectangle, amount, -3), RectangleInt.Inflate(original, amount, -3));
            var expected = rectangle;
            var actual = original;
            expected.Inflate(new Size(amount, -3));
            actual.Inflate(new SizeInt(amount, -3));
            Same(expected, actual);
            expected.Offset(new Point(amount, 2));
            actual.Offset(new PointInt(amount, 2));
            Same(expected, actual);
            expected.Offset(4, amount);
            actual.Offset(4, amount);
            Same(expected, actual);
            Same(rectangle, original);
        }
    }

    [Fact]
    public void FloatOffsetAndInflationMatchSystemDrawing()
    {
        foreach (var rectangle in FloatCases())
        foreach (float amount in new[] { -20.25f, 0, 7.5f, float.MaxValue, float.NaN })
        {
            RectangleFloat original = rectangle;
            Same(RectangleF.Inflate(rectangle, amount, -3.5f), RectangleFloat.Inflate(original, amount, -3.5f));
            var expected = rectangle;
            var actual = original;
            expected.Inflate(new SizeF(amount, -3.5f));
            actual.Inflate(new SizeFloat(amount, -3.5f));
            Same(expected, actual);
            expected.Offset(new PointF(amount, 2.25f));
            actual.Offset(new PointFloat(amount, 2.25f));
            Same(expected, actual);
            expected.Offset(4.5f, amount);
            actual.Offset(4.5f, amount);
            Same(expected, actual);
            Same(rectangle, original);
        }
    }

    [Theory]
    [InlineData(-2.5f, 3.5f, 4.5f, -5.5f)]
    [InlineData(-2.9f, 3.1f, 0.1f, -0.1f)]
    [InlineData(0, 0, 0, 0)]
    [InlineData(16777215, -16777215, 100000, 200000)]
    public void FloatToIntegerRoundingMatchesSystemDrawing(float x, float y, float width, float height)
    {
        var drawing = new RectangleF(x, y, width, height);
        var actual = new RectangleFloat(x, y, width, height);
        Same(Rectangle.Ceiling(drawing), RectangleInt.Ceiling(actual));
        Same(Rectangle.Round(drawing), RectangleInt.Round(actual));
        Same(Rectangle.Truncate(drawing), RectangleInt.Truncate(actual));
    }

    [Fact]
    public void ConversionsPreserveComponents()
    {
        foreach (var rectangle in IntegerCases())
        {
            RectangleInt integer = rectangle;
            Rectangle drawing = integer;
            Assert.Equal(rectangle, drawing);
            RectangleFloat floating = integer;
            Same((RectangleF)rectangle, floating);
        }
        foreach (var rectangle in FloatCases())
        {
            RectangleFloat floating = rectangle;
            RectangleF drawing = floating;
            Same(drawing, floating);
            var vector = new Vector4(rectangle.X, rectangle.Y, rectangle.Width, rectangle.Height);
            Same(rectangle, new RectangleFloat(vector));
            Same(rectangle, (RectangleFloat)vector);
            Assert.Equal(vector, floating.ToVector4());
            Assert.Equal(vector, (Vector4)floating);
        }
    }

    [Theory]
    [InlineData("en-US")]
    [InlineData("ru-RU")]
    public void StringRepresentationUsesCurrentCultureLikeSystemDrawing(string culture)
    {
        var previous = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo(culture);
            Assert.Equal(new Rectangle(-3, 4, 5, 6).ToString(), new RectangleInt(-3, 4, 5, 6).ToString());
            Assert.Equal(new RectangleF(-3.25f, 4.5f, 5.75f, 6.125f).ToString(),
                new RectangleFloat(-3.25f, 4.5f, 5.75f, 6.125f).ToString());
        }
        finally { CultureInfo.CurrentCulture = previous; }
    }

    [Fact]
    public void IntegerArithmeticRetainsImagingRoundingAndDoesNotMutateInputs()
    {
        var rectangle = new RectangleInt(10, 20, 3, 4);
        var point = new PointInt(-2, 5);
        Assert.Equal(new RectangleInt(8, 25, 3, 4), rectangle.Add(point));
        Assert.Equal(new RectangleInt(12, 15, 3, 4), rectangle.Sub(point));
        Assert.Equal(rectangle.Add(point), rectangle + point);
        Assert.Equal(rectangle.Sub(point), rectangle - point);
        Assert.Equal(12, rectangle.Area());
        Assert.Equal(12, RectangleInt.Area(rectangle.Size));
        Assert.Equal(new PointInt(10, 20), rectangle.GetPoint());
        Assert.Equal(new RectangleInt(10, 20, 4, 4), rectangle.ToBox());
        Assert.Equal(new RectangleInt(9, 19, 4, 6), rectangle.ToBox(.5f));
        Assert.Equal(new RectangleInt(10, 19, 4, 6), rectangle.Scale(.5f, .5f));
        Assert.Equal(new RectangleInt(10, 20, 4, 4), rectangle.Scale(kx: .5f));
        Assert.Equal(new RectangleInt(10, 19, 3, 6), rectangle.Scale(ky: .5f));
        Assert.Equal(new RectangleInt(9, 20, 5, 5), rectangle.Scale());
        Assert.Equal(rectangle, rectangle.Scale(0, 0));
        Assert.Equal(rectangle, rectangle.ToBox(0));
        Assert.Equal(new RectangleInt(10, 20, 3, 4), rectangle);

        var negativeOrigin = new RectangleInt(-2, -3, 3, 5);
        Assert.Equal(new RectangleInt(-2, -4, 4, 7), negativeOrigin.ToBox(.5f));
        Assert.Equal(new RectangleInt(-2, -4, 4, 7), negativeOrigin.Scale(.5f, .5f));
        Assert.Equal(new RectangleInt(-15000, -10000, 50000, 50000),
            new RectangleInt(-5000, -5000, 30000, 40000).Scale());
    }

    [Fact]
    public void FloatArithmeticPreservesFractionalCoordinates()
    {
        var rectangle = new RectangleFloat(10, 20, 3, 4);
        var point = new PointFloat(-2.25f, 5.5f);
        Assert.Equal(new RectangleFloat(7.75f, 25.5f, 3, 4), rectangle.Add(point));
        Assert.Equal(new RectangleFloat(12.25f, 14.5f, 3, 4), rectangle.Sub(point));
        Assert.Equal(rectangle.Add(point), rectangle + point);
        Assert.Equal(rectangle.Sub(point), rectangle - point);
        Assert.Equal(12, rectangle.Area());
        Assert.Equal(12, RectangleFloat.Area(rectangle.Size));
        Assert.Equal(new PointFloat(10, 20), rectangle.GetPoint());
        Assert.Equal(new RectangleFloat(9.5f, 20, 4, 4), rectangle.ToBox());
        Assert.Equal(new RectangleFloat(9.25f, 19, 4.5f, 6), rectangle.ToBox(.5f));
        Assert.Equal(new RectangleFloat(9.25f, 19, 4.5f, 6), rectangle.Scale(.5f, .5f));
        Assert.Equal(new RectangleFloat(9.25f, 20, 4.5f, 4), rectangle.Scale(kx: .5f));
        Assert.Equal(new RectangleFloat(10, 19, 3, 6), rectangle.Scale(ky: .5f));
        Assert.Equal(new RectangleFloat(9, 19.5f, 5, 5), rectangle.Scale());
        Assert.Equal(rectangle, rectangle.Scale(0, 0));
        Assert.Equal(rectangle, rectangle.ToBox(0));
        Assert.Equal(new RectangleFloat(10, 20, 3, 4), rectangle);
        NumericAssert.Close(5e20, new RectangleFloat(0, 0, 3e20f, 4e20f).Scale().Width);
    }

    [Fact]
    public void CornerOrderAndRoundTripMatchImagingConvention()
    {
        var integer = new RectangleInt(2, 3, 4, 5);
        Assert.Equal(new[] { new PointInt(2, 3), new PointInt(6, 3), new PointInt(6, 8), new PointInt(2, 8) },
            integer.ToPoints());
        Assert.Equal(integer, RectangleInt.FromPoints(integer.ToPoints()));
        var floating = new RectangleFloat(2.5f, 3.25f, 4.5f, 5.75f);
        Assert.Equal(new[] { new PointFloat(2.5f, 3.25f), new PointFloat(7, 3.25f),
            new PointFloat(7, 9), new PointFloat(2.5f, 9) }, floating.ToPoints());
        Assert.Equal(floating, RectangleFloat.FromPoints(floating.ToPoints()));
        Assert.Equal(RectangleInt.Empty, RectangleInt.FromPoints(RectangleInt.Empty.ToPoints()));
        var negative = new RectangleFloat(4, 5, -2, -3);
        Assert.Equal(negative, RectangleFloat.FromPoints(negative.ToPoints()));
        Assert.Throws<ArgumentNullException>(() => RectangleInt.FromPoints(null!));
        Assert.Throws<ArgumentNullException>(() => RectangleFloat.FromPoints(null!));
        foreach (int length in new[] { 0, 1, 3, 5 })
        {
            Assert.Throws<ArgumentException>(() => RectangleInt.FromPoints(new PointInt[length]));
            Assert.Throws<ArgumentException>(() => RectangleFloat.FromPoints(new PointFloat[length]));
        }
    }

    [Fact]
    public void ArrayOperationsReturnNewArraysAndLeaveInputsUnchanged()
    {
        var integers = new[] { new RectangleInt(10, 20, 3, 4), new RectangleInt(-2, -3, 3, 5) };
        var originals = integers.ToArray();
        Assert.Equal(new[] { new RectangleInt(11, 18, 3, 4), new RectangleInt(-1, -5, 3, 5) },
            RectangleInt.Add(integers, new PointInt(1, -2)));
        Assert.Equal(new[] { new RectangleInt(9, 22, 3, 4), new RectangleInt(-3, -1, 3, 5) },
            RectangleInt.Sub(integers, new PointInt(1, -2)));
        Assert.Equal(new[] { new RectangleInt(10, 20, 4, 4), new RectangleInt(-3, -3, 5, 5) },
            RectangleInt.ToBox(integers));
        Assert.Equal(new[] { new RectangleInt(9, 19, 4, 6), new RectangleInt(-2, -4, 4, 7) },
            RectangleInt.ToBox(.5f, integers));
        Assert.Equal(originals, integers);
        Assert.NotSame(integers, RectangleInt.Add(integers, default));

        var floats = new[] { new RectangleFloat(10, 20, 3, 4), new RectangleFloat(-2, -3, 3, 5) };
        var floatOriginals = floats.ToArray();
        Assert.Equal(new[] { new RectangleFloat(11, 18, 3, 4), new RectangleFloat(-1, -5, 3, 5) },
            RectangleFloat.Add(floats, new PointFloat(1, -2)));
        Assert.Equal(new[] { new RectangleFloat(9, 22, 3, 4), new RectangleFloat(-3, -1, 3, 5) },
            RectangleFloat.Sub(floats, new PointFloat(1, -2)));
        Assert.Equal(new[] { new RectangleFloat(9.5f, 20, 4, 4), new RectangleFloat(-3, -3, 5, 5) },
            RectangleFloat.ToBox(floats));
        Assert.Equal(new[] { new RectangleFloat(9.25f, 19, 4.5f, 6), new RectangleFloat(-2.75f, -4.25f, 4.5f, 7.5f) },
            RectangleFloat.ToBox(.5f, floats));
        Assert.Equal(floatOriginals, floats);
        Assert.NotSame(floats, RectangleFloat.Sub(floats, default));
        Assert.Empty(RectangleInt.ToBox(Array.Empty<RectangleInt>()));
        Assert.Empty(RectangleFloat.ToBox(.5f, Array.Empty<RectangleFloat>()));
    }

    [Fact]
    public void ArrayOperationsRejectNullInputs()
    {
        Assert.Throws<ArgumentNullException>(() => RectangleInt.Add(null!, default));
        Assert.Throws<ArgumentNullException>(() => RectangleInt.Sub(null!, default));
        Assert.Throws<ArgumentNullException>(() => RectangleInt.Min(null!));
        Assert.Throws<ArgumentNullException>(() => RectangleInt.Max(null!));
        Assert.Throws<ArgumentNullException>(() => RectangleInt.ToBox(null!));
        Assert.Throws<ArgumentNullException>(() => RectangleInt.ToBox(.5f, null!));
        Assert.Throws<ArgumentNullException>(() => RectangleFloat.Add(null!, default));
        Assert.Throws<ArgumentNullException>(() => RectangleFloat.Sub(null!, default));
        Assert.Throws<ArgumentNullException>(() => RectangleFloat.Min(null!));
        Assert.Throws<ArgumentNullException>(() => RectangleFloat.Max(null!));
        Assert.Throws<ArgumentNullException>(() => RectangleFloat.ToBox(null!));
        Assert.Throws<ArgumentNullException>(() => RectangleFloat.ToBox(.5f, null!));
    }

    [Fact]
    public void AreaSelectionHandlesEmptyInputsTiesAndOverflow()
    {
        var small = new RectangleInt(1, 2, 3, 4);
        var tie = new RectangleInt(9, 8, 6, 2);
        var large = new RectangleInt(1, 2, 100000, 100000);
        Assert.Equal(RectangleInt.Empty, RectangleInt.Min());
        Assert.Equal(RectangleInt.Empty, RectangleInt.Max());
        Assert.Equal(RectangleInt.Empty, RectangleInt.Min(RectangleInt.Empty));
        Assert.Equal(RectangleInt.Empty, RectangleInt.Max(RectangleInt.Empty));
        Assert.Equal(small, RectangleInt.Min(RectangleInt.Empty, small, tie, large));
        Assert.Equal(large, RectangleInt.Max(RectangleInt.Empty, small, large));
        Assert.Equal(small, RectangleInt.Max(small, tie));
        var maxArea = new RectangleInt(0, 0, int.MaxValue, 1);
        Assert.Equal(maxArea, RectangleInt.Min(maxArea));
        Assert.Equal(small, RectangleInt.Min(small, new RectangleInt(1, 1, 0, 0)));

        var tinyFloat = new RectangleFloat(2, 3, float.Epsilon, float.Epsilon);
        var hugeFloat = new RectangleFloat(2, 3, float.MaxValue, float.MaxValue);
        var smallerFloat = new RectangleFloat(4, 5, float.MaxValue, float.MaxValue / 2);
        Assert.Equal(RectangleFloat.Empty, RectangleFloat.Min());
        Assert.Equal(RectangleFloat.Empty, RectangleFloat.Max());
        Assert.Equal(RectangleFloat.Empty, RectangleFloat.Min(new RectangleFloat(3, 4, 0, 7)));
        Assert.Equal(RectangleFloat.Empty, RectangleFloat.Max(new RectangleFloat(3, 4, -1, 7)));
        Assert.Equal(tinyFloat, RectangleFloat.Min(hugeFloat, RectangleFloat.Empty, tinyFloat));
        Assert.Equal(hugeFloat, RectangleFloat.Max(smallerFloat, hugeFloat));
        Assert.Equal(smallerFloat, RectangleFloat.Min(hugeFloat, smallerFloat));
        Assert.Equal((RectangleFloat)small, RectangleFloat.Min(small, tie));
        Assert.Equal((RectangleFloat)small, RectangleFloat.Max(small, tie));
    }

    [Fact]
    public void AreaSelectionSkipsZeroAndNegativeDimensionsForBothTypes()
    {
        var empty = new[] { new RectangleInt(10, 20, 0, 5), new RectangleInt(10, 20, 5, 0),
            new RectangleInt(10, 20, -100, 100), new RectangleInt(10, 20, 100, -100),
            new RectangleInt(10, 20, -100, -100) };
        var valid = new RectangleInt(1, 2, 3, 4);
        foreach (var rectangle in empty)
        {
            Assert.Equal(RectangleInt.Empty, RectangleInt.Min(rectangle));
            Assert.Equal(RectangleInt.Empty, RectangleInt.Max(rectangle));
            Assert.Equal(valid, RectangleInt.Min(rectangle, valid));
            Assert.Equal(valid, RectangleInt.Max(rectangle, valid));

            RectangleFloat floating = rectangle;
            Assert.Equal(RectangleFloat.Empty, RectangleFloat.Min(floating));
            Assert.Equal(RectangleFloat.Empty, RectangleFloat.Max(floating));
            Assert.Equal((RectangleFloat)valid, RectangleFloat.Min(floating, valid));
            Assert.Equal((RectangleFloat)valid, RectangleFloat.Max(floating, valid));
        }
    }

    [Fact]
    public void IoUHandlesOverlapContainmentTouchingAndDegenerateRectangles()
    {
        var a = new RectangleInt(0, 0, 10, 10);
        var b = new RectangleInt(5, 5, 10, 10);
        NumericAssert.Close(1.0 / 7, a.IoU(b));
        NumericAssert.Close(1.0 / 7, RectangleInt.IoU(b, a));
        Assert.Equal(1, a.IoU(a));
        Assert.Equal(.25f, a.IoU(new RectangleInt(2, 2, 5, 5)));
        Assert.Equal(0, a.IoU(new RectangleInt(10, 0, 10, 10)));
        Assert.Equal(0, a.IoU(new RectangleInt(50, 50, 10, 10)));
        Assert.Equal(0, a.IoU(RectangleInt.Empty));
        Assert.Equal(0, RectangleInt.Empty.IoU(RectangleInt.Empty));
        Assert.Equal(0, a.IoU(new RectangleInt(5, 5, -2, 3)));

        var af = new RectangleFloat(.5f, .25f, 10, 10);
        var bf = new RectangleFloat(5.5f, 5.25f, 10, 10);
        NumericAssert.Close(1.0 / 7, af.IoU(bf));
        NumericAssert.Close(1.0 / 7, RectangleFloat.IoU(bf, af));
        Assert.Equal(1, af.IoU(af));
        Assert.Equal(0, af.IoU(new RectangleFloat(10.5f, .25f, 10, 10)));
        Assert.Equal(0, af.IoU(RectangleFloat.Empty));
        Assert.Equal(0, RectangleFloat.Empty.IoU(RectangleFloat.Empty));
        Assert.Equal(0, af.IoU(new RectangleFloat(5, 5, 2, -3)));
    }

    [Fact]
    public void IoUAvoidsCoordinateAreaOverflowAndUnderflow()
    {
        Assert.Equal(1, new RectangleInt(0, 0, 100000, 100000).IoU(new RectangleInt(0, 0, 100000, 100000)));
        var edge = new RectangleInt(int.MaxValue - 2, int.MinValue, 10, 10);
        Assert.Equal(1, edge.IoU(edge));
        var other = new RectangleInt(int.MaxValue - 7, int.MinValue, 10, 10);
        NumericAssert.Close(1.0 / 3, edge.IoU(other));
        Assert.Equal(0, edge.IoU(new RectangleInt(int.MinValue, int.MinValue, 10, 10)));
        var enormous = new RectangleFloat(0, 0, float.MaxValue, float.MaxValue);
        Assert.Equal(1, enormous.IoU(enormous));
        var tiny = new RectangleFloat(0, 0, float.Epsilon, float.Epsilon);
        Assert.Equal(1, tiny.IoU(tiny));
    }

    [Fact]
    public void ClampNormalizesOnlyTheFirstRectangleAndRetainsDisjointOrigins()
    {
        var bounds = new RectangleInt(0, 0, 10, 10);
        var reversed = new RectangleInt(12, 8, -8, -12);
        Assert.Equal(new RectangleInt(4, 0, 6, 8), reversed.Clamp(bounds));
        Assert.Equal(new RectangleInt(4, 0, 6, 8), RectangleInt.Clamp(reversed, bounds));
        Assert.Equal(new RectangleInt(20, 2, 0, 3), new RectangleInt(20, 2, 5, 3).Clamp(bounds));
        Assert.Equal(new RectangleInt(0, 2, 0, 3), new RectangleInt(-20, 2, 5, 3).Clamp(bounds));
        Assert.Equal(new RectangleInt(12, 8, -8, -12), reversed);
        Assert.Equal(new RectangleInt(5, 5, 0, 0), bounds.Clamp(new RectangleInt(5, 5, -3, -3)));
        Assert.Equal(new RectangleInt(-10, 0, 10, 1),
            new RectangleInt(0, 0, int.MinValue, 1).Clamp(new RectangleInt(-10, 0, 20, 1)));
        Assert.Equal(new RectangleInt(int.MaxValue, 0, 5, 1),
            new RectangleInt(int.MaxValue, 0, 10, 1).Clamp(new RectangleInt(int.MaxValue, 0, 5, 1)));

        var floatBounds = new RectangleFloat(0, 0, 10, 10);
        var floatReversed = new RectangleFloat(12.5f, 8.25f, -8, -12);
        Assert.Equal(new RectangleFloat(4.5f, 0, 5.5f, 8.25f), floatReversed.Clamp(floatBounds));
        Assert.Equal(new RectangleFloat(4.5f, 0, 5.5f, 8.25f), RectangleFloat.Clamp(floatReversed, floatBounds));
        Assert.Equal(new RectangleFloat(20, 2, 0, 3), new RectangleFloat(20, 2, 5, 3).Clamp(floatBounds));
        Assert.Equal(new RectangleFloat(5, 5, 0, 0), floatBounds.Clamp(new RectangleFloat(5, 5, -3, -3)));
        Assert.Equal(new RectangleFloat(12.5f, 8.25f, -8, -12), floatReversed);
    }
}
