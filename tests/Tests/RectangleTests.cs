using System.Drawing;
using System.Globalization;
using System.Numerics;
using UMapx.Core;
using Xunit;
using ReferenceRectangles = UMapx.Tests.Reference.Rectangles;

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
        SameBits(expected.X, actual.X);
        SameBits(expected.Y, actual.Y);
        SameBits(expected.Width, actual.Width);
        SameBits(expected.Height, actual.Height);
    }

    private static Rectangle[] IntegerCases(int seed = 731)
    {
        var cases = new List<Rectangle>
        {
            Rectangle.Empty, new(3, 4, 0, 0), new(0, 0, 10, 20), new(10, 20, 5, 5),
            new(10, 0, 4, 3), new(-4, -3, 10, 20), new(2, 3, -4, 5), new(2, 3, 4, -5),
            new(2, 3, -4, -5), new(5, 5, 0, 4), new(5, 5, 4, 0),
            new(int.MaxValue, int.MinValue, 4, -5), new(int.MinValue, int.MaxValue, -4, 5),
            new(-1, -2, int.MaxValue, int.MinValue)
        };
        var random = new Random(seed);
        for (int i = 0; i < 60; i++)
            cases.Add(new Rectangle(random.Next(-50, 50), random.Next(-50, 50),
                random.Next(-10, 70), random.Next(-10, 70)));
        for (int i = 0; i < 64; i++)
            cases.Add(new Rectangle((int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1),
                (int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1),
                (int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1),
                (int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1)));
        return cases.ToArray();
    }

    private static RectangleF[] FloatCases(int seed = 731)
    {
        var cases = IntegerCases(seed).Select(r => (RectangleF)r).ToList();
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
        var random = new Random(seed);
        for (int i = 0; i < 64; i++)
        {
            float Next() => BitConverter.Int32BitsToSingle((int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1));
            cases.Add(new RectangleF(Next(), Next(), Next(), Next()));
        }
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
            Assert.Equal(expected.GetHashCode(), actual.GetHashCode());

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
            Assert.Equal(expected.GetHashCode(), actual.GetHashCode());

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
            RectangleInt actual = Convert(a);
            RectangleInt other = Convert(b);
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
            RectangleFloat actual = Convert(a);
            RectangleFloat other = Convert(b);
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
            RectangleInt original = Convert(rectangle);
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
            RectangleFloat original = Convert(rectangle);
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
            RectangleInt integer = Convert(rectangle);
            RectangleFloat floating = integer;
            Same((RectangleF)rectangle, floating);
        }
        foreach (var rectangle in FloatCases())
        {
            RectangleFloat floating = Convert(rectangle);
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

    private static RectangleInt Convert(Rectangle r) => new RectangleInt(r.X, r.Y, r.Width, r.Height);
    private static RectangleFloat Convert(RectangleF r) => new RectangleFloat(r.X, r.Y, r.Width, r.Height);
    private static PointInt Convert(Point p) => new PointInt(p.X, p.Y);
    private static PointFloat Convert(PointF p) => new PointFloat(p.X, p.Y);

    private static void SameBits(float expected, float actual) =>
        Assert.Equal(BitConverter.SingleToInt32Bits(expected), BitConverter.SingleToInt32Bits(actual));

    private static void Same(Point expected, PointInt actual)
    {
        Assert.Equal(expected.X, actual.X);
        Assert.Equal(expected.Y, actual.Y);
    }

    private static void Same(PointF expected, PointFloat actual)
    {
        SameBits(expected.X, actual.X);
        SameBits(expected.Y, actual.Y);
    }

    private static void SameArray(Rectangle[] expected, RectangleInt[] actual)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++) Same(expected[i], actual[i]);
    }

    private static void SameArray(RectangleF[] expected, RectangleFloat[] actual)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++) Same(expected[i], actual[i]);
    }

    private static void SameResult<TReference, TActual>(
        Func<TReference> reference, Func<TActual> actual, Action<TReference, TActual> compare)
    {
        TReference expectedValue = default!;
        TActual actualValue = default!;
        Exception? expectedException = null;
        Exception? actualException = null;
        try { expectedValue = reference(); } catch (Exception exception) { expectedException = exception; }
        try { actualValue = actual(); } catch (Exception exception) { actualException = exception; }
        Assert.Equal(expectedException?.GetType(), actualException?.GetType());
        if (expectedException is ArgumentException argument)
            Assert.Equal(argument.Message, actualException!.Message);
        if (expectedException == null) compare(expectedValue, actualValue);
    }

    private static readonly float[] Factors = { -2, -1, -0.5f, -0.0f, 0, 0.5f, 1, 2.75f,
        float.Epsilon, float.MaxValue, float.MinValue, float.NaN,
        float.PositiveInfinity, float.NegativeInfinity };

    [Fact]
    public void RoundingExtremeFloatComponentsMatchesSystemDrawing()
    {
        foreach (var rectangle in FloatCases())
        {
            var actual = Convert(rectangle);
            Same(Rectangle.Ceiling(rectangle), RectangleInt.Ceiling(actual));
            Same(Rectangle.Round(rectangle), RectangleInt.Round(actual));
            Same(Rectangle.Truncate(rectangle), RectangleInt.Truncate(actual));
        }
    }

    [Theory]
    [InlineData(731)]
    [InlineData(1999)]
    [InlineData(42)]
    public void IntegerArithmeticMatchesReferenceIncludingOverflowAndSpecialValues(int seed)
    {
        var cases = IntegerCases(seed);
        foreach (var rectangle in cases)
        {
            var actual = Convert(rectangle);
            Same(ReferenceRectangles.ToBox(rectangle), actual.ToBox());
            Same(ReferenceRectangles.Scale(rectangle), actual.Scale());
            Same(ReferenceRectangles.GetPoint(rectangle), actual.GetPoint());
            Assert.Equal(ReferenceRectangles.Area(rectangle), actual.Area());
            Assert.Equal(ReferenceRectangles.Area(rectangle.Size), RectangleInt.Area(actual.Size));
            var corners = ReferenceRectangles.ToPoints(rectangle);
            var actualCorners = actual.ToPoints();
            Assert.Equal(corners.Length, actualCorners.Length);
            for (int i = 0; i < corners.Length; i++) Same(corners[i], actualCorners[i]);
            Same(ReferenceRectangles.FromPoints(corners), RectangleInt.FromPoints(actualCorners));
            foreach (float factor in Factors)
            {
                Same(ReferenceRectangles.ToBox(rectangle, factor), actual.ToBox(factor));
                Same(ReferenceRectangles.Scale(rectangle, kx: factor), actual.Scale(kx: factor));
                Same(ReferenceRectangles.Scale(rectangle, ky: factor), actual.Scale(ky: factor));
                Same(ReferenceRectangles.Scale(rectangle, factor, -factor), actual.Scale(factor, -factor));
            }
            foreach (var other in cases)
            {
                Same(ReferenceRectangles.Add(rectangle, other.Location), actual.Add(Convert(other.Location)));
                Same(ReferenceRectangles.Sub(rectangle, other.Location), actual.Sub(Convert(other.Location)));
                SameBits(ReferenceRectangles.IoU(rectangle, other), actual.IoU(Convert(other)));
                Same(ReferenceRectangles.Clamp(rectangle, other), actual.Clamp(Convert(other)));
            }
            Same(rectangle, actual);
        }
    }

    [Theory]
    [InlineData(731)]
    [InlineData(1999)]
    [InlineData(42)]
    public void IntegerArrayArithmeticMatchesReferenceWithAgreedEmptinessRule(int seed)
    {
        var pool = IntegerCases(seed);
        var arrays = new List<Rectangle[]?> { null, Array.Empty<Rectangle>(), pool,
            new[] { new Rectangle(3, 4, 0, 0), new Rectangle(1, 2, -3, -4) },
            new[] { new Rectangle(1, 2, 3, 4), new Rectangle(5, 6, 4, 3) } };
        for (int i = 0; i < pool.Length; i++)
        {
            arrays.Add(new[] { pool[i] });
            arrays.Add(new[] { pool[i], pool[(i * 7 + 1) % pool.Length], pool[(i * 13 + 2) % pool.Length] });
        }
        for (int i = 0; i < arrays.Count; i++)
        {
            var rectangles = arrays[i];
            var actual = rectangles?.Select(Convert).ToArray();
            var offset = pool[i % pool.Length].Location;
            SameResult(() => ReferenceRectangles.Add(rectangles!, offset),
                () => RectangleInt.Add(actual!, Convert(offset)), SameArray);
            SameResult(() => ReferenceRectangles.Sub(rectangles!, offset),
                () => RectangleInt.Sub(actual!, Convert(offset)), SameArray);
            SameResult(() => ReferenceRectangles.Max(rectangles!), () => RectangleInt.Max(actual!), Same);
            SameResult(() => ReferenceRectangles.Min(rectangles!), () => RectangleInt.Min(actual!), Same);
            SameResult(() => ReferenceRectangles.ToBox(rectangles!), () => RectangleInt.ToBox(actual!), SameArray);
            foreach (float factor in Factors)
                SameResult(() => ReferenceRectangles.ToBox(factor, rectangles!),
                    () => RectangleInt.ToBox(factor, actual!), SameArray);
            if (rectangles != null)
            {
                SameArray(rectangles, actual!);
                Assert.NotSame(actual, RectangleInt.Add(actual!, Convert(Point.Empty)));
                Assert.NotSame(actual, RectangleInt.Sub(actual!, Convert(Point.Empty)));
                Assert.NotSame(actual, RectangleInt.ToBox(actual!));
                Assert.NotSame(actual, RectangleInt.ToBox(0.5f, actual!));
            }
        }
    }

    [Fact]
    public void IntegerFromPointsMatchesReferenceForUnorderedAndInvalidInputs()
    {
        var pool = IntegerCases().Select(r => r.Location).ToArray();
        var arrays = new List<Point[]?> { null, Array.Empty<Point>(), pool };
        for (int length = 1; length <= 5; length++)
        for (int i = 0; i < pool.Length; i++)
            arrays.Add(Enumerable.Range(0, length).Select(n => pool[(i + 11 * n) % pool.Length]).ToArray());
        foreach (var points in arrays)
            SameResult(() => ReferenceRectangles.FromPoints(points!),
                () => RectangleInt.FromPoints(points?.Select(Convert).ToArray()!), Same);
    }

    [Theory]
    [InlineData(731)]
    [InlineData(1999)]
    [InlineData(42)]
    public void FloatArithmeticMatchesReferenceIncludingOverflowAndSpecialValues(int seed)
    {
        var cases = FloatCases(seed);
        foreach (var rectangle in cases)
        {
            var actual = Convert(rectangle);
            Same(ReferenceRectangles.ToBox(rectangle), actual.ToBox());
            Same(ReferenceRectangles.Scale(rectangle), actual.Scale());
            Same(ReferenceRectangles.GetPoint(rectangle), actual.GetPoint());
            SameBits(ReferenceRectangles.Area(rectangle), actual.Area());
            SameBits(ReferenceRectangles.Area(rectangle.Size), RectangleFloat.Area(actual.Size));
            var corners = ReferenceRectangles.ToPoints(rectangle);
            var actualCorners = actual.ToPoints();
            Assert.Equal(corners.Length, actualCorners.Length);
            for (int i = 0; i < corners.Length; i++) Same(corners[i], actualCorners[i]);
            Same(ReferenceRectangles.FromPoints(corners), RectangleFloat.FromPoints(actualCorners));
            foreach (float factor in Factors)
            {
                Same(ReferenceRectangles.ToBox(rectangle, factor), actual.ToBox(factor));
                Same(ReferenceRectangles.Scale(rectangle, kx: factor), actual.Scale(kx: factor));
                Same(ReferenceRectangles.Scale(rectangle, ky: factor), actual.Scale(ky: factor));
                Same(ReferenceRectangles.Scale(rectangle, factor, -factor), actual.Scale(factor, -factor));
            }
            foreach (var other in cases)
            {
                Same(ReferenceRectangles.Add(rectangle, other.Location), actual.Add(Convert(other.Location)));
                Same(ReferenceRectangles.Sub(rectangle, other.Location), actual.Sub(Convert(other.Location)));
                SameBits(ReferenceRectangles.IoU(rectangle, other), actual.IoU(Convert(other)));
                Same(ReferenceRectangles.Clamp(rectangle, other), actual.Clamp(Convert(other)));
            }
            Same(rectangle, actual);
        }
    }

    [Theory]
    [InlineData(731)]
    [InlineData(1999)]
    [InlineData(42)]
    public void FloatArrayArithmeticMatchesReferenceWithAgreedEmptinessRule(int seed)
    {
        var pool = FloatCases(seed);
        var arrays = new List<RectangleF[]?> { null, Array.Empty<RectangleF>(), pool,
            new[] { new RectangleF(3, 4, 0, 0), new RectangleF(1, 2, -3, -4) },
            new[] { new RectangleF(1, 2, 3, 4), new RectangleF(5, 6, 4, 3) } };
        for (int i = 0; i < pool.Length; i++)
        {
            arrays.Add(new[] { pool[i] });
            arrays.Add(new[] { pool[i], pool[(i * 7 + 1) % pool.Length], pool[(i * 13 + 2) % pool.Length] });
        }
        for (int i = 0; i < arrays.Count; i++)
        {
            var rectangles = arrays[i];
            var actual = rectangles?.Select(Convert).ToArray();
            var offset = pool[i % pool.Length].Location;
            SameResult(() => ReferenceRectangles.Add(rectangles!, offset),
                () => RectangleFloat.Add(actual!, Convert(offset)), SameArray);
            SameResult(() => ReferenceRectangles.Sub(rectangles!, offset),
                () => RectangleFloat.Sub(actual!, Convert(offset)), SameArray);
            SameResult(() => ReferenceRectangles.Max(rectangles!), () => RectangleFloat.Max(actual!), Same);
            SameResult(() => ReferenceRectangles.Min(rectangles!), () => RectangleFloat.Min(actual!), Same);
            SameResult(() => ReferenceRectangles.ToBox(rectangles!), () => RectangleFloat.ToBox(actual!), SameArray);
            foreach (float factor in Factors)
                SameResult(() => ReferenceRectangles.ToBox(factor, rectangles!),
                    () => RectangleFloat.ToBox(factor, actual!), SameArray);
            if (rectangles != null)
            {
                SameArray(rectangles, actual!);
                Assert.NotSame(actual, RectangleFloat.Add(actual!, Convert(PointF.Empty)));
                Assert.NotSame(actual, RectangleFloat.Sub(actual!, Convert(PointF.Empty)));
                Assert.NotSame(actual, RectangleFloat.ToBox(actual!));
                Assert.NotSame(actual, RectangleFloat.ToBox(0.5f, actual!));
            }
        }
    }

    [Fact]
    public void FloatFromPointsMatchesReferenceForUnorderedAndInvalidInputs()
    {
        var pool = FloatCases().Select(r => r.Location).ToArray();
        var arrays = new List<PointF[]?> { null, Array.Empty<PointF>(), pool };
        for (int length = 1; length <= 5; length++)
        for (int i = 0; i < pool.Length; i++)
            arrays.Add(Enumerable.Range(0, length).Select(n => pool[(i + 11 * n) % pool.Length]).ToArray());
        foreach (var points in arrays)
            SameResult(() => ReferenceRectangles.FromPoints(points!),
                () => RectangleFloat.FromPoints(points?.Select(Convert).ToArray()!), Same);
    }
}
