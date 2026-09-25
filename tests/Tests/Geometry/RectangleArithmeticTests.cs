using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Geometry")]
public class RectangleArithmeticTests
{
    [Fact]
    public void IntegerTranslationAndBoxOperationsPreserveRoundingAndInput()
    {
        var rectangle = new RectangleInt(10, 20, 3, 4);
        var offset = new PointInt(-2, 5);
        Assert.Equal(new RectangleInt(8, 25, 3, 4), rectangle.Add(offset));
        Assert.Equal(new RectangleInt(12, 15, 3, 4), rectangle.Sub(offset));
        Assert.Equal(new RectangleInt(10, 20, 4, 4), rectangle.ToBox());
        Assert.Equal(new RectangleInt(9, 19, 4, 6), rectangle.ToBox(0.5f));
        Assert.Equal(new RectangleInt(10, 19, 4, 6), rectangle.Scale(0.5f, 0.5f));
        Assert.Equal(new RectangleInt(9, 20, 5, 5), rectangle.Scale());
        Assert.Equal(rectangle, rectangle.Scale(0, 0));
        Assert.Equal(new RectangleInt(10, 20, 3, 4), rectangle);
    }

    [Fact]
    public void FloatBoxOperationsPreserveTheirIntegerTruncation()
    {
        var rectangle = new RectangleFloat(10.25f, 20.5f, 3, 4);
        var offset = new PointFloat(-2.5f, 5.25f);
        Assert.Equal(new RectangleFloat(7.75f, 25.75f, 3, 4), rectangle.Add(offset));
        Assert.Equal(new RectangleFloat(12.75f, 15.25f, 3, 4), rectangle.Sub(offset));
        Assert.Equal(new RectangleFloat(9.75f, 20.5f, 4, 4), rectangle.ToBox());
        Assert.Equal(new RectangleFloat(9, 19, 4, 6), rectangle.ToBox(0.5f));
        Assert.Equal(new RectangleFloat(10.25f, 19.5f, 4, 6), rectangle.Scale(0.5f, 0.5f));
        Assert.Equal(new RectangleFloat(9.25f, 20, 5, 5), rectangle.Scale());
        Assert.Equal(rectangle, rectangle.Scale(0, 0));
        Assert.Equal(new RectangleFloat(10.25f, 20.5f, 3, 4), rectangle);
    }

    [Fact]
    public void IntegerArrayOperationsPreserveOrderAndInput()
    {
        var rectangles = new[]
        {
            new RectangleInt(10, 20, 3, 4),
            new RectangleInt(-10, -20, 6, 2)
        };
        var original = (RectangleInt[])rectangles.Clone();
        var offset = new PointInt(-2, 5);
        var added = RectangleInt.Add(rectangles, offset);
        var subtracted = RectangleInt.Sub(rectangles, offset);
        var boxes = RectangleInt.ToBox(rectangles);
        var enlarged = RectangleInt.ToBox(0.5f, rectangles);
        Assert.Equal(new[] { new RectangleInt(8, 25, 3, 4), new RectangleInt(-12, -15, 6, 2) }, added);
        Assert.Equal(new[] { new RectangleInt(12, 15, 3, 4), new RectangleInt(-8, -25, 6, 2) }, subtracted);
        Assert.Equal(new[] { new RectangleInt(10, 20, 4, 4), new RectangleInt(-10, -22, 6, 6) }, boxes);
        Assert.Equal(new[] { new RectangleInt(9, 19, 4, 6), new RectangleInt(-11, -20, 9, 3) }, enlarged);
        Assert.NotSame(rectangles, added);
        Assert.NotSame(rectangles, subtracted);
        Assert.NotSame(rectangles, boxes);
        Assert.NotSame(rectangles, enlarged);
        Assert.Equal(original, rectangles);
    }

    [Fact]
    public void FloatArrayOperationsPreserveOrderAndInput()
    {
        var rectangles = new[]
        {
            new RectangleFloat(10.25f, 20.5f, 3, 4),
            new RectangleFloat(-10.25f, -20.5f, 6, 2)
        };
        var original = (RectangleFloat[])rectangles.Clone();
        var offset = new PointFloat(-2.5f, 5.25f);
        var added = RectangleFloat.Add(rectangles, offset);
        var subtracted = RectangleFloat.Sub(rectangles, offset);
        var boxes = RectangleFloat.ToBox(rectangles);
        var enlarged = RectangleFloat.ToBox(0.5f, rectangles);
        Assert.Equal(new[] { new RectangleFloat(7.75f, 25.75f, 3, 4), new RectangleFloat(-12.75f, -15.25f, 6, 2) }, added);
        Assert.Equal(new[] { new RectangleFloat(12.75f, 15.25f, 3, 4), new RectangleFloat(-7.75f, -25.75f, 6, 2) }, subtracted);
        Assert.Equal(new[] { new RectangleFloat(9.75f, 20.5f, 4, 4), new RectangleFloat(-10.25f, -22.5f, 6, 6) }, boxes);
        Assert.Equal(new[] { new RectangleFloat(9, 19, 4, 6), new RectangleFloat(-11, -21, 9, 3) }, enlarged);
        Assert.NotSame(rectangles, added);
        Assert.NotSame(rectangles, subtracted);
        Assert.NotSame(rectangles, boxes);
        Assert.NotSame(rectangles, enlarged);
        Assert.Equal(original, rectangles);
    }

    [Fact]
    public void CornersAreClockwiseAndFromPointsUsesOnlyFirstAndThirdPoints()
    {
        var rectangle = new RectangleInt(10, 20, 3, 4);
        Assert.Equal(new PointInt(10, 20), rectangle.GetPoint());
        Assert.Equal(new[] { new PointInt(10, 20), new PointInt(13, 20), new PointInt(13, 24), new PointInt(10, 24) }, rectangle.ToPoints());
        Assert.Equal(rectangle, RectangleInt.FromPoints(new[] { new PointInt(10, 20), new PointInt(-99, 99), new PointInt(13, 24), new PointInt(99, -99) }));
        var floating = new RectangleFloat(10.25f, 20.5f, 3.5f, 4.25f);
        Assert.Equal(new PointFloat(10.25f, 20.5f), floating.GetPoint());
        Assert.Equal(new[] { new PointFloat(10.25f, 20.5f), new PointFloat(13.75f, 20.5f), new PointFloat(13.75f, 24.75f), new PointFloat(10.25f, 24.75f) }, floating.ToPoints());
        Assert.Equal(floating, RectangleFloat.FromPoints(new[] { new PointFloat(10.25f, 20.5f), new PointFloat(-99, 99), new PointFloat(13.75f, 24.75f), new PointFloat(99, -99) }));
    }

    [Theory]
    [InlineData(0)]
    [InlineData(3)]
    [InlineData(5)]
    public void FromPointsRejectsArraysWithoutExactlyFourPoints(int length)
    {
        Assert.Throws<ArgumentException>(() => RectangleInt.FromPoints(new PointInt[length]));
        Assert.Throws<ArgumentException>(() => RectangleFloat.FromPoints(new PointFloat[length]));
    }

    [Fact]
    public void AreaPreservesSignAndIntegerOverflow()
    {
        Assert.Equal(-12, new RectangleInt(10, 20, -3, 4).Area());
        Assert.Equal(-12, RectangleInt.Area(new SizeInt(-3, 4)));
        Assert.Equal(-2, new RectangleInt(0, 0, int.MaxValue, 2).Area());
        Assert.Equal(-2, RectangleInt.Area(new SizeInt(int.MaxValue, 2)));
        Assert.Equal(-14, new RectangleFloat(10, 20, -3.5f, 4).Area());
        Assert.Equal(-14, RectangleFloat.Area(new SizeFloat(-3.5f, 4)));
    }

    [Fact]
    public void MinAndMaxSkipEmptyRectanglesAndKeepTheFirstTie()
    {
        var empty = new RectangleInt(99, 98, -4, -5);
        var small = new RectangleInt(1, 2, 1, 2);
        var large = new RectangleInt(3, 4, 3, 4);
        var tie = new RectangleInt(5, 6, 4, 3);
        Assert.Equal(small, RectangleInt.Min(empty, large, small));
        Assert.Equal(large, RectangleInt.Max(empty, small, large, tie));
        Assert.Equal(large, RectangleInt.Min(large, tie));
        Assert.Equal(empty, RectangleInt.Max(empty));
        Assert.Equal(RectangleInt.Empty, RectangleInt.Min(empty));
        Assert.Equal(RectangleInt.Empty, RectangleInt.Max());
        Assert.Equal(RectangleInt.Empty, RectangleInt.Min());
        var floatingEmpty = new RectangleFloat(99, 98, -4, -5);
        var floatingSmall = new RectangleFloat(1, 2, 0.5f, 2);
        var floatingLarge = new RectangleFloat(3, 4, 3.5f, 4);
        var floatingTie = new RectangleFloat(5, 6, 4, 3.5f);
        Assert.Equal(floatingSmall, RectangleFloat.Min(floatingEmpty, floatingLarge, floatingSmall));
        Assert.Equal(floatingLarge, RectangleFloat.Max(floatingEmpty, floatingSmall, floatingLarge, floatingTie));
        Assert.Equal(floatingLarge, RectangleFloat.Min(floatingLarge, floatingTie));
        Assert.Equal(floatingEmpty, RectangleFloat.Max(floatingEmpty));
        Assert.Equal(RectangleFloat.Empty, RectangleFloat.Min(floatingEmpty));
        Assert.Equal(RectangleFloat.Empty, RectangleFloat.Max());
        Assert.Equal(RectangleFloat.Empty, RectangleFloat.Min());
    }

    [Theory]
    [InlineData(0, 1)]
    [InlineData(2, 1.0 / 3)]
    [InlineData(4, 0)]
    [InlineData(5, 0)]
    public void IoUHandlesCoincidentOverlappingTouchingAndDisjointRectangles(int x, double expected)
    {
        var rectangle = new RectangleInt(0, 0, 4, 4);
        var other = new RectangleInt(x, 0, 4, 4);
        NumericAssert.Close(expected, rectangle.IoU(other), 1e-7, 0);
        NumericAssert.Close(expected, other.IoU(rectangle), 1e-7, 0);
        var floating = new RectangleFloat(0.25f, 0.5f, 4, 4);
        var floatingOther = new RectangleFloat(x + 0.25f, 0.5f, 4, 4);
        NumericAssert.Close(expected, floating.IoU(floatingOther), 1e-7, 0);
        NumericAssert.Close(expected, floatingOther.IoU(floating), 1e-7, 0);
    }

    [Theory]
    [InlineData(12, 8, -8, -12, 4, 0, 6, 8)]
    [InlineData(20, 2, 5, 3, 20, 2, 0, 3)]
    [InlineData(2, 20, 3, 5, 2, 20, 3, 0)]
    [InlineData(2, 3, 4, 5, 2, 3, 4, 5)]
    public void ClampNormalizesNegativeDimensionsAndRetainsDisjointLocation(int x, int y, int width, int height, int expectedX, int expectedY, int expectedWidth, int expectedHeight)
    {
        var rectangle = new RectangleInt(x, y, width, height);
        Assert.Equal(new RectangleInt(expectedX, expectedY, expectedWidth, expectedHeight), rectangle.Clamp(new RectangleInt(0, 0, 10, 10)));
        Assert.Equal(new RectangleInt(x, y, width, height), rectangle);
        var floating = new RectangleFloat(x + 0.25f, y + 0.5f, width, height);
        Assert.Equal(new RectangleFloat(expectedX + 0.25f, expectedY + 0.5f, expectedWidth, expectedHeight), floating.Clamp(new RectangleFloat(0.25f, 0.5f, 10, 10)));
        Assert.Equal(new RectangleFloat(x + 0.25f, y + 0.5f, width, height), floating);
    }
}
