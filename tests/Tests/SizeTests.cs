using System.Drawing;
using System.Globalization;
using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class SizeTests
{
    private static readonly int[] Integers = { int.MinValue, int.MaxValue, -65536, -1, 0, 1, 2, 65536 };
    private static readonly float[] Floats = { float.MinValue, float.MaxValue, float.NaN,
        float.NegativeInfinity, float.PositiveInfinity, -float.Epsilon, float.Epsilon,
        -0.0f, 0, -3.5f, 2.5f, 16777216f };

    private static void Same(Size expected, SizeInt actual)
    {
        Assert.Equal(expected.Width, actual.Width);
        Assert.Equal(expected.Height, actual.Height);
    }

    private static void Same(float expected, float actual) =>
        Assert.Equal(BitConverter.SingleToInt32Bits(expected), BitConverter.SingleToInt32Bits(actual));

    private static void Same(SizeF expected, SizeFloat actual)
    {
        Same(expected.Width, actual.Width);
        Same(expected.Height, actual.Height);
    }

    private static void SameArithmetic(SizeF expected, SizeFloat actual)
    {
        if (float.IsNaN(expected.Width)) Assert.True(float.IsNaN(actual.Width));
        else Same(expected.Width, actual.Width);
        if (float.IsNaN(expected.Height)) Assert.True(float.IsNaN(actual.Height));
        else Same(expected.Height, actual.Height);
    }

    private static void SameResult<TExpected, TActual>(
        Func<TExpected> expected, Func<TActual> actual, Action<TExpected, TActual> compare)
    {
        TExpected expectedResult = default!;
        TActual actualResult = default!;
        Exception? expectedError = null, actualError = null;
        try { expectedResult = expected(); } catch (Exception error) { expectedError = error; }
        try { actualResult = actual(); } catch (Exception error) { actualError = error; }
        Assert.Equal(expectedError?.GetType(), actualError?.GetType());
        if (expectedError == null) compare(expectedResult, actualResult);
    }

    [Fact]
    public void IntegerConstructionPropertiesEqualityAndCloningMatchSystemDrawing()
    {
        Same(Size.Empty, SizeInt.Empty);
        foreach (int width in Integers)
        foreach (int height in Integers)
        {
            var expected = new Size(width, height);
            var actual = new SizeInt(width, height);
            Same(expected, actual);
            Same(new Size(new Point(width, height)), new SizeInt(new PointInt(width, height)));
            Assert.Equal(expected.IsEmpty, actual.IsEmpty);
            Assert.Equal(expected.GetHashCode(), actual.GetHashCode());
            Assert.Equal(expected.ToString(), actual.ToString());
            Assert.Equal(expected == Size.Empty, actual == SizeInt.Empty);
            Assert.Equal(expected != Size.Empty, actual != SizeInt.Empty);
            Assert.True(((IEquatable<SizeInt>)actual).Equals(actual));
            Assert.True(actual.Equals((object)actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals(expected));
            Same(expected, actual.Clone());
            Same(expected, (SizeInt)((ICloneable)actual).Clone());
            SizeFloat floating = actual;
            Same((SizeF)expected, floating);
            Assert.Equal(new PointInt(width, height), (PointInt)actual);
            var changed = actual;
            changed.Width = height;
            changed.Height = width;
            Same(new Size(height, width), changed);
            Same(expected, actual);
        }
    }

    [Fact]
    public void FloatConstructionPropertiesEqualityAndConversionsMatchSystemDrawing()
    {
        Same(SizeF.Empty, SizeFloat.Empty);
        foreach (float width in Floats)
        foreach (float height in Floats)
        {
            var expected = new SizeF(width, height);
            var actual = new SizeFloat(width, height);
            Same(expected, actual);
            Same(new SizeF(expected), new SizeFloat(actual));
            Same(new SizeF(new PointF(width, height)), new SizeFloat(new PointFloat(width, height)));
            var vector = new Vector2(width, height);
            Same(expected, new SizeFloat(vector));
            Same(expected, (SizeFloat)vector);
            Same(width, actual.ToVector2().X);
            Same(height, actual.ToVector2().Y);
            Same(width, ((Vector2)actual).X);
            Same(height, ((Vector2)actual).Y);
            Same(expected, actual.Clone());
            Same(expected, (SizeFloat)((ICloneable)actual).Clone());
            Assert.Equal(expected.IsEmpty, actual.IsEmpty);
            Assert.Equal(expected.GetHashCode(), actual.GetHashCode());
            Assert.Equal(expected.ToString(), actual.ToString());
            Assert.Equal(expected == SizeF.Empty, actual == SizeFloat.Empty);
            Assert.Equal(expected != SizeF.Empty, actual != SizeFloat.Empty);
            Assert.Equal(expected.Equals(expected), ((IEquatable<SizeFloat>)actual).Equals(actual));
            Assert.Equal(expected.Equals((object)expected), actual.Equals((object)actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals(expected));
            var point = actual.ToPointF();
            var cast = (PointFloat)actual;
            Same(expected.ToPointF().X, point.X);
            Same(expected.ToPointF().Y, point.Y);
            Same(point.X, cast.X);
            Same(point.Y, cast.Y);
            Same(expected.ToSize(), actual.ToSize());
            Same(Size.Ceiling(expected), SizeInt.Ceiling(actual));
            Same(Size.Round(expected), SizeInt.Round(actual));
            Same(Size.Truncate(expected), SizeInt.Truncate(actual));
            var changed = actual;
            changed.Width = height;
            changed.Height = width;
            Same(new SizeF(height, width), changed);
            Same(expected, actual);
        }
    }

    [Fact]
    public void IntegerAdditionAndSubtractionMatchSystemDrawing()
    {
        foreach (int width in Integers)
        foreach (int height in Integers)
        foreach (int secondWidth in Integers)
        foreach (int secondHeight in Integers)
        {
            var expected = new Size(width, height);
            var actual = new SizeInt(width, height);
            var second = new Size(secondWidth, secondHeight);
            var secondActual = new SizeInt(secondWidth, secondHeight);
            Same(Size.Add(expected, second), SizeInt.Add(actual, secondActual));
            Same(Size.Subtract(expected, second), SizeInt.Subtract(actual, secondActual));
            Same(expected + second, actual + secondActual);
            Same(expected - second, actual - secondActual);
            Assert.Equal(expected.Equals(second), actual.Equals(secondActual));
            Assert.Equal(expected == second, actual == secondActual);
            Assert.Equal(expected != second, actual != secondActual);
            Same(new Size(width, height), actual);
            Same(new Size(secondWidth, secondHeight), secondActual);
        }
    }

    [Fact]
    public void IntegerScalingMatchesSystemDrawingIncludingDivisionExceptions()
    {
        foreach (int width in Integers)
        foreach (int height in Integers)
        {
            var expected = new Size(width, height);
            var actual = new SizeInt(width, height);
            foreach (int factor in Integers)
            {
                Same(expected * factor, actual * factor);
                Same(factor * expected, factor * actual);
                SameResult(() => expected / factor, () => actual / factor, Same);
            }
            foreach (float floating in Floats)
            {
                Same(expected * floating, actual * floating);
                Same(floating * expected, floating * actual);
                Same(expected / floating, actual / floating);
            }
            Same(new Size(width, height), actual);
        }
    }

    [Theory]
    [InlineData(731)]
    [InlineData(42)]
    public void FloatArithmeticAndRoundingMatchSystemDrawing(int seed)
    {
        var random = new Random(seed);
        var sizes = new List<SizeF>();
        foreach (float width in Floats)
        foreach (float height in Floats)
            sizes.Add(new SizeF(width, height));
        for (int i = 0; i < 256; i++)
        {
            float Next() => BitConverter.Int32BitsToSingle(
                (int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1));
            sizes.Add(new SizeF(Next(), Next()));
        }
        foreach (var expected in sizes)
        {
            var actual = new SizeFloat(expected.Width, expected.Height);
            foreach (var second in sizes)
            {
                var secondActual = new SizeFloat(second.Width, second.Height);
                SameArithmetic(SizeF.Add(expected, second), SizeFloat.Add(actual, secondActual));
                SameArithmetic(SizeF.Subtract(expected, second), SizeFloat.Subtract(actual, secondActual));
                SameArithmetic(expected + second, actual + secondActual);
                SameArithmetic(expected - second, actual - secondActual);
                SameArithmetic(expected * second.Width, actual * second.Width);
                SameArithmetic(second.Width * expected, second.Width * actual);
                SameArithmetic(expected / second.Width, actual / second.Width);
                Assert.Equal(expected.Equals(second), actual.Equals(secondActual));
                Assert.Equal(expected == second, actual == secondActual);
                Assert.Equal(expected != second, actual != secondActual);
            }
            Same(Size.Ceiling(expected), SizeInt.Ceiling(actual));
            Same(Size.Round(expected), SizeInt.Round(actual));
            Same(Size.Truncate(expected), SizeInt.Truncate(actual));
            Same(expected, actual);
        }
    }

    [Theory]
    [InlineData("en-US")]
    [InlineData("ru-RU")]
    public void FormattingMatchesSystemDrawing(string culture)
    {
        var previous = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo(culture);
            Assert.Equal(new Size(-3, 7).ToString(), new SizeInt(-3, 7).ToString());
            Assert.Equal(new SizeF(-3.25f, 7.5f).ToString(), new SizeFloat(-3.25f, 7.5f).ToString());
        }
        finally { CultureInfo.CurrentCulture = previous; }
    }
}
