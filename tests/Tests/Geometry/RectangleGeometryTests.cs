using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Geometry")]
public class RectangleGeometryTests
{
    [Fact]
    public void GeometryMatchesEnumeratedGridPointsAndCells()
    {
        var cases = (
            from x in Enumerable.Range(-2, 5)
            from y in Enumerable.Range(-2, 5)
            from width in Enumerable.Range(1, 3)
            from height in Enumerable.Range(1, 3)
            select new
            {
                Bounds = (x, y, width, height),
                Cells = Grid(x, y, width, height),
                Vertices = Grid(x, y, width + 1, height + 1)
            }

        ).ToArray();
        foreach (var a in cases)
            foreach (var b in cases)
            {
                var integer = Integer(a.Bounds);
                var otherInteger = Integer(b.Bounds);
                var floating = Floating(a.Bounds, 0.5f);
                var otherFloat = Floating(b.Bounds, 0.5f);
                var intersection = Bounds(a.Vertices.Intersect(b.Vertices));
                var union = Bounds(a.Vertices.Union(b.Vertices));
                bool contains = b.Vertices.IsSubsetOf(a.Vertices);
                bool overlaps = a.Cells.Overlaps(b.Cells);
                bool containsPoint = a.Cells.Contains((b.Bounds.x, b.Bounds.y));
                bool equal = a.Bounds == b.Bounds;
                Assert.Equal(equal, integer == otherInteger);
                Assert.Equal(!equal, integer != otherInteger);
                Assert.Equal(equal, integer.Equals(otherInteger));
                Assert.Equal(equal, floating == otherFloat);
                Assert.Equal(!equal, floating != otherFloat);
                Assert.Equal(equal, floating.Equals(otherFloat));
                Assert.Equal(contains, integer.Contains(otherInteger));
                Assert.Equal(contains, floating.Contains(otherFloat));
                Assert.Equal(overlaps, integer.IntersectsWith(otherInteger));
                Assert.Equal(overlaps, floating.IntersectsWith(otherFloat));
                Assert.Equal(containsPoint, integer.Contains(b.Bounds.x, b.Bounds.y));
                Assert.Equal(containsPoint, integer.Contains(new PointInt(b.Bounds.x, b.Bounds.y)));
                Assert.Equal(containsPoint, floating.Contains(b.Bounds.x * 0.5f, b.Bounds.y * 0.5f));
                Assert.Equal(containsPoint, floating.Contains(new PointFloat(b.Bounds.x * 0.5f, b.Bounds.y * 0.5f)));
                Assert.Equal(Integer(intersection), RectangleInt.Intersect(integer, otherInteger));
                Assert.Equal(Floating(intersection, 0.5f), RectangleFloat.Intersect(floating, otherFloat));
                Assert.Equal(Integer(union), RectangleInt.Union(integer, otherInteger));
                Assert.Equal(Floating(union, 0.5f), RectangleFloat.Union(floating, otherFloat));
                Assert.Equal(Integer(a.Bounds), integer);
                Assert.Equal(Floating(a.Bounds, 0.5f), floating);
                integer.Intersect(otherInteger);
                floating.Intersect(otherFloat);
                Assert.Equal(Integer(intersection), integer);
                Assert.Equal(Floating(intersection, 0.5f), floating);
            }
    }

    private static HashSet<(int X, int Y)> Grid(int x, int y, int width, int height) => (
        from column in Enumerable.Range(x, width) from row in Enumerable.Range(y, height) select (column, row)).ToHashSet();
    private static (int X, int Y, int Width, int Height) Bounds(IEnumerable<(int X, int Y)> points)
    {
        var values = points.ToArray();
        if (values.Length == 0)
            return (0, 0, 0, 0);
        int left = values.Min(p => p.X), top = values.Min(p => p.Y);
        return (left, top, values.Max(p => p.X) - left, values.Max(p => p.Y) - top);
    }

    private static RectangleInt Integer((int X, int Y, int Width, int Height) r) => new(r.X, r.Y, r.Width, r.Height);
    private static RectangleFloat Floating((int X, int Y, int Width, int Height) r, float scale = 1) => new(r.X * scale, r.Y * scale, r.Width * scale, r.Height * scale);
    [Theory]
    [InlineData(5, 2, 0, 4, true, true, 5, 2, 0, 4, 10, 10)]
    [InlineData(10, 2, 0, 4, true, false, 10, 2, 0, 4, 10, 10)]
    [InlineData(5, 2, -2, 4, true, true, 0, 0, 0, 0, 10, 10)]
    [InlineData(2, 3, -4, 5, true, false, 0, 0, 0, 0, 10, 10)]
    [InlineData(2, 5, 4, -2, true, true, 0, 0, 0, 0, 10, 10)]
    [InlineData(20, 2, 0, 4, false, false, 0, 0, 0, 0, 20, 10)]
    [InlineData(0, 0, 0, 0, true, false, 0, 0, 0, 0, 10, 10)]
    [InlineData(3, 4, 0, 0, true, true, 3, 4, 0, 0, 10, 10)]
    [InlineData(10, 10, 4, 4, false, false, 10, 10, 0, 0, 14, 14)]
    public void EmptyNegativeAndTouchingBoundsPreserveSpecifiedGeometry(int x, int y, int width, int height, bool contains, bool overlaps, int ix, int iy, int iw, int ih, int uw, int uh)
    {
        var a = new RectangleInt(0, 0, 10, 10);
        var b = new RectangleInt(x, y, width, height);
        Assert.Equal(contains, a.Contains(b));
        Assert.Equal(overlaps, a.IntersectsWith(b));
        Assert.Equal(overlaps, b.IntersectsWith(a));
        var expected = new RectangleInt(ix, iy, iw, ih);
        Assert.Equal(expected, RectangleInt.Intersect(a, b));
        Assert.Equal(expected, RectangleInt.Intersect(b, a));
        Assert.Equal(new RectangleInt(0, 0, uw, uh), RectangleInt.Union(a, b));
        Assert.Equal(new RectangleInt(0, 0, uw, uh), RectangleInt.Union(b, a));
        a.Intersect(b);
        Assert.Equal(expected, a);
        var af = new RectangleFloat(0, 0, 10, 10);
        var bf = new RectangleFloat(x, y, width, height);
        Assert.Equal(contains, af.Contains(bf));
        Assert.Equal(overlaps, af.IntersectsWith(bf));
        Assert.Equal(overlaps, bf.IntersectsWith(af));
        var expectedFloat = new RectangleFloat(ix, iy, iw, ih);
        Assert.Equal(expectedFloat, RectangleFloat.Intersect(af, bf));
        Assert.Equal(expectedFloat, RectangleFloat.Intersect(bf, af));
        Assert.Equal(new RectangleFloat(0, 0, uw, uh), RectangleFloat.Union(af, bf));
        Assert.Equal(new RectangleFloat(0, 0, uw, uh), RectangleFloat.Union(bf, af));
        af.Intersect(bf);
        Assert.Equal(expectedFloat, af);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void IntegerGeometryUsesWrappedBounds(bool transpose)
    {
        RectangleInt Create(int x, int width) => transpose ? new RectangleInt(0, x, 10, width) : new RectangleInt(x, 0, width, 10);
        var normal = new RectangleInt(0, 0, 10, 10);
        var overflow = Create(int.MaxValue - 1, 4);
        Assert.Equal(int.MinValue + 2, transpose ? overflow.Bottom : overflow.Right);
        Assert.False(overflow.Contains(0, 0));
        Assert.False(overflow.Contains(normal));
        Assert.True(normal.Contains(overflow));
        Assert.False(overflow.IntersectsWith(normal));
        Assert.False(normal.IntersectsWith(overflow));
        Assert.Equal(RectangleInt.Empty, RectangleInt.Intersect(overflow, normal));
        Assert.Equal(normal, RectangleInt.Union(overflow, normal));
        var underflow = Create(int.MinValue, -1);
        Assert.Equal(int.MaxValue, transpose ? underflow.Bottom : underflow.Right);
        Assert.True(underflow.IsEmpty);
        Assert.True(underflow.Contains(0, 0));
        Assert.True(underflow.Contains(normal));
        Assert.False(normal.Contains(underflow));
        Assert.True(underflow.IntersectsWith(normal));
        Assert.True(normal.IntersectsWith(underflow));
        Assert.Equal(normal, RectangleInt.Intersect(underflow, normal));
        Assert.Equal(underflow, RectangleInt.Union(underflow, normal));
        overflow.Intersect(normal);
        underflow.Intersect(normal);
        Assert.Equal(RectangleInt.Empty, overflow);
        Assert.Equal(normal, underflow);
    }

    [Theory]
    [InlineData(0, float.NaN, 0, float.NaN, 10)]
    [InlineData(1, 0, float.NaN, 10, float.NaN)]
    [InlineData(2, 0, 0, float.NaN, 10)]
    [InlineData(3, 0, 0, 10, float.NaN)]
    public void NaNComponentsPreventContainmentAndPropagateThroughUnion(int component, float x, float y, float width, float height)
    {
        float[] values =
        {
            2,
            3,
            4,
            5
        };
        values[component] = float.NaN;
        var invalid = new RectangleFloat(values[0], values[1], values[2], values[3]);
        var normal = new RectangleFloat(0, 0, 10, 10);
        Assert.False(invalid.Contains(2, 3));
        Assert.False(invalid.Contains(new PointFloat(2, 3)));
        Assert.False(invalid.Contains(normal));
        Assert.False(normal.Contains(invalid));
        Assert.False(invalid.IntersectsWith(normal));
        Assert.False(normal.IntersectsWith(invalid));
        Assert.Equal(RectangleFloat.Empty, RectangleFloat.Intersect(invalid, normal));
        Assert.Equal(RectangleFloat.Empty, RectangleFloat.Intersect(normal, invalid));
        Same((x, y, width, height), RectangleFloat.Union(invalid, normal));
        Same((x, y, width, height), RectangleFloat.Union(normal, invalid));
        var copy = invalid;
        Assert.False(invalid == copy);
        Assert.True(invalid != copy);
        Assert.False(invalid.Equals(copy));
        invalid.Intersect(normal);
        Assert.Equal(RectangleFloat.Empty, invalid);
    }

    [Theory]
    [InlineData(2, float.PositiveInfinity, true, 2, 3, 8, 4, 0, float.PositiveInfinity)]
    [InlineData(float.PositiveInfinity, 4, false, 0, 0, 0, 0, 0, float.PositiveInfinity)]
    [InlineData(float.NegativeInfinity, float.PositiveInfinity, false, 0, 0, 0, 0, float.NegativeInfinity, float.NaN)]
    [InlineData(float.MaxValue, float.MaxValue, false, 0, 0, 0, 0, 0, float.PositiveInfinity)]
    [InlineData(float.NegativeInfinity, 4, false, 0, 0, 0, 0, float.NegativeInfinity, float.PositiveInfinity)]
    public void InfiniteBoundsHaveSpecifiedIntersectionsAndUnions(float x, float width, bool overlaps, float ix, float iy, float iw, float ih, float ux, float uw)
    {
        var a = new RectangleFloat(0, 0, 10, 10);
        var b = new RectangleFloat(x, 3, width, 4);
        Assert.False(a.Contains(b));
        Assert.False(b.Contains(a));
        Assert.Equal(overlaps, a.IntersectsWith(b));
        Assert.Equal(overlaps, b.IntersectsWith(a));
        Same((ix, iy, iw, ih), RectangleFloat.Intersect(a, b));
        Same((ix, iy, iw, ih), RectangleFloat.Intersect(b, a));
        Same((ux, 0, uw, 10), RectangleFloat.Union(a, b));
        Same((ux, 0, uw, 10), RectangleFloat.Union(b, a));
        a.Intersect(b);
        Same((ix, iy, iw, ih), a);
    }

    [Fact]
    public void SignedZeroBoundsPreserveExtremaSigns()
    {
        var negative = new RectangleFloat(-0.0f, -0.0f, -0.0f, -0.0f);
        var union = RectangleFloat.Union(RectangleFloat.Empty, negative);
        var intersection = RectangleFloat.Intersect(RectangleFloat.Empty, negative);
        Assert.Equal(int.MinValue, BitConverter.SingleToInt32Bits(union.X));
        Assert.Equal(int.MinValue, BitConverter.SingleToInt32Bits(union.Y));
        Assert.Equal(0, BitConverter.SingleToInt32Bits(union.Width));
        Assert.Equal(0, BitConverter.SingleToInt32Bits(union.Height));
        Assert.Equal(0, BitConverter.SingleToInt32Bits(intersection.X));
        Assert.Equal(0, BitConverter.SingleToInt32Bits(intersection.Y));
        Assert.Equal(int.MinValue, BitConverter.SingleToInt32Bits(intersection.Width));
        Assert.Equal(int.MinValue, BitConverter.SingleToInt32Bits(intersection.Height));
    }

    private static void Same((float X, float Y, float Width, float Height) expected, RectangleFloat actual)
    {
        Assert.Equal(expected.X, actual.X);
        Assert.Equal(expected.Y, actual.Y);
        Assert.Equal(expected.Width, actual.Width);
        Assert.Equal(expected.Height, actual.Height);
    }
}
