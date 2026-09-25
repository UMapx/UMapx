using System.Drawing;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class RangeTests
{
    private static void Same(RangeInt expected, RangeInt actual)
    {
        Assert.Equal(expected.Min, actual.Min);
        Assert.Equal(expected.Max, actual.Max);
    }

    private static void Same(float expected, float actual) =>
        Assert.Equal(BitConverter.SingleToInt32Bits(expected), BitConverter.SingleToInt32Bits(actual));

    private static void Same(RangeFloat expected, RangeFloat actual)
    {
        Same(expected.Min, actual.Min);
        Same(expected.Max, actual.Max);
    }

    [Fact]
    public void IntegerOperationsMatchClosedSetsIncludingTouchingAndReversedBounds()
    {
        for (int min = -3; min <= 3; min++)
        for (int max = -3; max <= 3; max++)
        {
            var a = new RangeInt(min, max);
            var first = Enumerable.Range(-3, 7).Where(x => min <= x && x <= max).ToArray();
            Assert.Equal(max - min, a.Length);
            foreach (int x in Enumerable.Range(-4, 9))
            {
                Assert.Equal(first.Contains(x), a.Contains(x));
                Assert.Equal(first.Contains(x), a.IsOnRange(x));
            }
            for (int otherMin = -3; otherMin <= 3; otherMin++)
            for (int otherMax = -3; otherMax <= 3; otherMax++)
            {
                var b = new RangeInt(otherMin, otherMax);
                var second = Enumerable.Range(-3, 7).Where(x => otherMin <= x && x <= otherMax).ToArray();
                var common = first.Intersect(second).ToArray();
                var combined = first.Concat(second).ToArray();
                var expectedIntersection = common.Length == 0 ? RangeInt.Empty :
                    new RangeInt(common.Min(), common.Max());
                var expectedUnion = combined.Length == 0 ? RangeInt.Empty :
                    new RangeInt(combined.Min(), combined.Max());
                Assert.Equal(common.Length > 0, a.IntersectsWith(b));
                Assert.Equal(first.Length > 0 && second.Length > 0 && second.All(first.Contains), a.Contains(b));
                Same(expectedIntersection, RangeInt.Intersect(a, b));
                Same(expectedUnion, RangeInt.Union(a, b));
                var changed = a;
                changed.Intersect(b);
                Same(expectedIntersection, changed);
                Same(new RangeInt(min, max), a);
            }
        }
    }

    [Fact]
    public void FloatOperationsMatchClosedSetsAtFractionalBounds()
    {
        var values = Enumerable.Range(-6, 13).Select(x => x / 2f).ToArray();
        foreach (float min in values)
        foreach (float max in values)
        {
            var a = new RangeFloat(min, max);
            var first = values.Where(x => min <= x && x <= max).ToArray();
            Same(max - min, a.Length);
            foreach (float x in values)
            {
                Assert.Equal(first.Contains(x), a.Contains(x));
                Assert.Equal(first.Contains(x), a.IsOnRange(x));
            }
            foreach (float otherMin in values)
            foreach (float otherMax in values)
            {
                var b = new RangeFloat(otherMin, otherMax);
                var second = values.Where(x => otherMin <= x && x <= otherMax).ToArray();
                var common = first.Intersect(second).ToArray();
                var combined = first.Concat(second).ToArray();
                var expectedIntersection = common.Length == 0 ? RangeFloat.Empty :
                    new RangeFloat(common.Min(), common.Max());
                var expectedUnion = combined.Length == 0 ? RangeFloat.Empty :
                    new RangeFloat(combined.Min(), combined.Max());
                Assert.Equal(common.Length > 0, a.IntersectsWith(b));
                Assert.Equal(first.Length > 0 && second.Length > 0 && second.All(first.Contains), a.Contains(b));
                Same(expectedIntersection, RangeFloat.Intersect(a, b));
                Same(expectedUnion, RangeFloat.Union(a, b));
                var changed = a;
                changed.Intersect(b);
                Same(expectedIntersection, changed);
                Same(new RangeFloat(min, max), a);
            }
        }
    }

    [Fact]
    public void EmptyRangeContainsNoValuesAndDefaultRangeContainsZero()
    {
        Same(new RangeInt(0, -1), RangeInt.Empty);
        Same(new RangeFloat(0, -1), RangeFloat.Empty);
        Assert.False(RangeInt.Empty.Contains(0));
        Assert.False(RangeFloat.Empty.Contains(0));
        Assert.True(default(RangeInt).Contains(0));
        Assert.True(default(RangeFloat).Contains(0));
        Assert.Null(typeof(RangeInt).GetProperty("IsEmpty"));
        Assert.Null(typeof(RangeFloat).GetProperty("IsEmpty"));
    }

    [Fact]
    public void FloatOperationsHandleNaNInfinityAndSignedZero()
    {
        var finite = new RangeFloat(-1, 1);
        foreach (var invalid in new[] { new RangeFloat(float.NaN, 1), new RangeFloat(0, float.NaN) })
        {
            Assert.False(invalid.Contains(0));
            Assert.False(invalid.Contains(finite));
            Assert.False(finite.Contains(invalid));
            Assert.False(invalid.IntersectsWith(finite));
            Assert.False(finite.IntersectsWith(invalid));
            Same(RangeFloat.Empty, RangeFloat.Intersect(invalid, finite));
            Same(RangeFloat.Empty, RangeFloat.Intersect(finite, invalid));
            Same(finite, RangeFloat.Union(invalid, finite));
            Same(finite, RangeFloat.Union(finite, invalid));
        }
        var all = new RangeFloat(float.NegativeInfinity, float.PositiveInfinity);
        Assert.True(all.Contains(finite));
        Assert.True(all.Contains(float.NegativeInfinity));
        Assert.True(all.Contains(float.PositiveInfinity));
        Assert.False(all.Contains(float.NaN));
        Same(finite, RangeFloat.Intersect(all, finite));
        Same(all, RangeFloat.Union(all, finite));
        var zero = new RangeFloat(-0.0f, 0.0f);
        Assert.True(zero.Contains(0));
        Same(zero, zero.Clone());
        Same(zero, (RangeFloat)((ICloneable)zero).Clone());
        Same(float.PositiveInfinity, all.Length);
    }

    [Fact]
    public void OffsetsAndInflationPreserveArithmeticAndDoNotNormalizeBounds()
    {
        var integer = new RangeInt(3, 1);
        integer.Offset(-2);
        Same(new RangeInt(1, -1), integer);
        integer.Inflate(2);
        Same(new RangeInt(-1, 1), integer);
        Same(new RangeInt(1, -1), RangeInt.Inflate(integer, -2));
        Same(new RangeInt(-1, 1), integer);
        var extremes = new RangeInt(int.MinValue, int.MaxValue);
        Assert.Equal(-1, extremes.Length);
        extremes.Offset(1);
        Same(new RangeInt(int.MinValue + 1, int.MinValue), extremes);
        extremes = new RangeInt(int.MinValue, int.MaxValue);
        extremes.Inflate(1);
        Same(new RangeInt(int.MaxValue, int.MinValue), extremes);

        var floating = new RangeFloat(3.5f, 1.5f);
        floating.Offset(-2.25f);
        Same(new RangeFloat(1.25f, -0.75f), floating);
        floating.Inflate(2.5f);
        Same(new RangeFloat(-1.25f, 1.75f), floating);
        Same(new RangeFloat(1.75f, -1.25f), RangeFloat.Inflate(floating, -3));
        Same(new RangeFloat(-1.25f, 1.75f), floating);
    }

    [Theory]
    [InlineData(-2.5f, 3.5f)]
    [InlineData(3.75f, -1.25f)]
    [InlineData(float.NaN, float.PositiveInfinity)]
    [InlineData(float.MinValue, float.MaxValue)]
    public void RoundingBoundsUsesSystemDrawingRounding(float min, float max)
    {
        var point = new PointF(min, max);
        var range = new RangeFloat(min, max);
        var ceiling = Point.Ceiling(point);
        var round = Point.Round(point);
        var truncate = Point.Truncate(point);
        Same(new RangeInt(ceiling.X, ceiling.Y), RangeInt.Ceiling(range));
        Same(new RangeInt(round.X, round.Y), RangeInt.Round(range));
        Same(new RangeInt(truncate.X, truncate.Y), RangeInt.Truncate(range));
    }

    [Fact]
    public void EqualityCloningAndConversionPreserveBounds()
    {
        var integer = new RangeInt(int.MinValue, int.MaxValue);
        Assert.True(((IEquatable<RangeInt>)integer).Equals(integer.Clone()));
        Assert.True(integer.Equals((object)integer.Clone()));
        Assert.False(integer.Equals(null));
        Same(integer, (RangeInt)((ICloneable)integer).Clone());
        RangeFloat floating = integer;
        Same((float)integer.Min, floating.Min);
        Same((float)integer.Max, floating.Max);
        Assert.True(((IEquatable<RangeFloat>)floating).Equals(floating.Clone()));
        Assert.True(floating.Equals((object)floating.Clone()));
        Assert.False(floating.Equals(null));
        Same(floating, (RangeFloat)((ICloneable)floating).Clone());
        var clone = floating.Clone();
        clone.Min = 10;
        clone.Max = 20;
        Same((float)integer.Min, floating.Min);
        Same((float)integer.Max, floating.Max);
        var nan = new RangeFloat(float.NaN, 1);
        Assert.False(((IEquatable<RangeFloat>)nan).Equals(nan));
        Assert.False(nan.Equals((object)nan));
    }
}
