using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class PointArithmeticTests
{
    [Fact]
    public void IntegerTranslationReturnsNewArraysAndPreservesOrderAndInput()
    {
        var points = new[] { new PointInt(-3, 5), new PointInt(2, -2) };
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
        var points = new[] { new PointFloat(-3.25f, 5.5f), new PointFloat(2.5f, -2.25f) };
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
        var points = new[] { new PointInt(3, 2), center };
        var rotated = PointInt.Rotate(points, center, 90);
        Assert.Equal(new[] { new PointInt(0, 3), center }, rotated);
        Assert.Equal(new PointInt(0, 3), points[0].Rotate(center, 90));
        Assert.Equal(new PointInt(1, 1), new PointInt(2, 0).Rotate(PointInt.Empty, 60));
        Assert.NotSame(points, rotated);
        Assert.Equal(new[] { new PointInt(3, 2), center }, points);

        var floatingCenter = new PointFloat(1.25f, 1.5f);
        var floatingPoints = new[] { new PointFloat(3.25f, 2.5f), floatingCenter };
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
        var points = new[] { new PointInt(-3, 5), new PointInt(2, -2) };
        Assert.Equal(new RectangleInt(-3, -2, 5, 7), PointInt.GetRectangle(points));
        Assert.Equal(new PointInt(0, 1), PointInt.GetMeanPoint(points));
        Assert.Equal(new PointInt(2, 5), points[0].GetSupportedPoint(points[1]));
        var floating = new[] { new PointFloat(-3.25f, 5.5f), new PointFloat(2.5f, -2.25f) };
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
        Assert.Equal(new RectangleFloat(float.MaxValue, float.MaxValue, float.NegativeInfinity, float.NegativeInfinity),
            PointFloat.GetRectangle(Array.Empty<PointFloat>()));
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
        Assert.Equal(new RectangleInt(int.MinValue, int.MinValue, -1, -1),
            PointInt.GetRectangle(new[] { point, new PointInt(int.MinValue, int.MaxValue) }));
    }
}
