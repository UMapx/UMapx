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

    private static void Same((int Width, int Height) expected, SizeInt actual)
    {
        Assert.Equal(expected.Width, actual.Width);
        Assert.Equal(expected.Height, actual.Height);
    }

    private static void Same(float expected, float actual) =>
        Assert.Equal(BitConverter.SingleToInt32Bits(expected), BitConverter.SingleToInt32Bits(actual));

    private static void Same((float Width, float Height) expected, SizeFloat actual)
    {
        Same(expected.Width, actual.Width);
        Same(expected.Height, actual.Height);
    }

    private static void SameArithmetic((float Width, float Height) expected, SizeFloat actual)
    {
        if (float.IsNaN(expected.Width)) Assert.True(float.IsNaN(actual.Width));
        else Same(expected.Width, actual.Width);
        if (float.IsNaN(expected.Height)) Assert.True(float.IsNaN(actual.Height));
        else Same(expected.Height, actual.Height);
    }

    [Fact]
    public void IntegerConstructionPropertiesEqualityAndCloningPreserveComponents()
    {
        Same((0, 0), SizeInt.Empty);
        foreach (int width in Integers)
        foreach (int height in Integers)
        {
            var expected = (Width: width, Height: height);
            var actual = new SizeInt(width, height);
            Same(expected, actual);
            Same(expected, new SizeInt(new PointInt(width, height)));
            Assert.Equal(width == 0 && height == 0, actual.IsEmpty);
            Assert.Equal(width == 0 && height == 0, actual == SizeInt.Empty);
            Assert.Equal(width != 0 || height != 0, actual != SizeInt.Empty);
            Assert.True(((IEquatable<SizeInt>)actual).Equals(actual));
            Assert.True(actual.Equals((object)actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals(expected));
            Same(expected, actual.Clone());
            Same(expected, (SizeInt)((ICloneable)actual).Clone());
            SizeFloat floating = actual;
            Same(((float)width, (float)height), floating);
            Assert.Equal(new PointInt(width, height), (PointInt)actual);
            var changed = actual;
            changed.Width = height;
            changed.Height = width;
            Same((height, width), changed);
            Same(expected, actual);
        }
    }

    [Fact]
    public void FloatConstructionPropertiesEqualityAndConversionsPreserveComponents()
    {
        Same((0f, 0f), SizeFloat.Empty);
        foreach (float width in Floats)
        foreach (float height in Floats)
        {
            var expected = (Width: width, Height: height);
            var actual = new SizeFloat(width, height);
            Same(expected, actual);
            Same(expected, new SizeFloat(actual));
            Same(expected, new SizeFloat(new PointFloat(width, height)));
            var vector = new Vector2(width, height);
            Same(expected, new SizeFloat(vector));
            Same(expected, (SizeFloat)vector);
            Same(width, actual.ToVector2().X);
            Same(height, actual.ToVector2().Y);
            Same(width, ((Vector2)actual).X);
            Same(height, ((Vector2)actual).Y);
            Same(expected, actual.Clone());
            Same(expected, (SizeFloat)((ICloneable)actual).Clone());
            Assert.Equal(width == 0 && height == 0, actual.IsEmpty);
            Assert.Equal(width == 0 && height == 0, actual == SizeFloat.Empty);
            Assert.Equal(width != 0 || height != 0, actual != SizeFloat.Empty);
            bool reflexive = !float.IsNaN(width) && !float.IsNaN(height);
            Assert.Equal(reflexive, ((IEquatable<SizeFloat>)actual).Equals(actual));
            Assert.Equal(reflexive, actual.Equals((object)actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals(expected));
            var point = actual.ToPointF();
            var cast = (PointFloat)actual;
            Same(width, point.X);
            Same(height, point.Y);
            Same(point.X, cast.X);
            Same(point.Y, cast.Y);
            Same((unchecked((int)width), unchecked((int)height)), actual.ToSize());
            Same((unchecked((int)Math.Ceiling(width)), unchecked((int)Math.Ceiling(height))), SizeInt.Ceiling(actual));
            Same((unchecked((int)Math.Round(width)), unchecked((int)Math.Round(height))), SizeInt.Round(actual));
            Same((unchecked((int)width), unchecked((int)height)), SizeInt.Truncate(actual));
            var changed = actual;
            changed.Width = height;
            changed.Height = width;
            Same((height, width), changed);
            Same(expected, actual);
        }
    }

    [Fact]
    public void IntegerAdditionAndSubtractionWrapOverflowInEachComponent()
    {
        foreach (int width in Integers)
        foreach (int height in Integers)
        foreach (int secondWidth in Integers)
        foreach (int secondHeight in Integers)
        {
            var actual = new SizeInt(width, height);
            var secondActual = new SizeInt(secondWidth, secondHeight);
            var sum = (unchecked(width + secondWidth), unchecked(height + secondHeight));
            var difference = (unchecked(width - secondWidth), unchecked(height - secondHeight));
            Same(sum, SizeInt.Add(actual, secondActual));
            Same(difference, SizeInt.Subtract(actual, secondActual));
            Same(sum, actual + secondActual);
            Same(difference, actual - secondActual);
            bool equal = width == secondWidth && height == secondHeight;
            Assert.Equal(equal, actual.Equals(secondActual));
            Assert.Equal(equal, actual == secondActual);
            Assert.Equal(!equal, actual != secondActual);
            Same((width, height), actual);
            Same((secondWidth, secondHeight), secondActual);
        }
    }

    [Fact]
    public void IntegerScalingUsesScalarArithmeticIncludingDivisionExceptions()
    {
        foreach (int width in Integers)
        foreach (int height in Integers)
        {
            var actual = new SizeInt(width, height);
            foreach (int factor in Integers)
            {
                var product = (unchecked(width * factor), unchecked(height * factor));
                Same(product, actual * factor);
                Same(product, factor * actual);
                if (factor == 0)
                    Assert.Throws<DivideByZeroException>(() => actual / factor);
                else if (factor == -1 && (width == int.MinValue || height == int.MinValue))
                    Assert.Throws<OverflowException>(() => actual / factor);
                else
                    Same((width / factor, height / factor), actual / factor);
            }
            foreach (float floating in Floats)
            {
                SameArithmetic((width * floating, height * floating), actual * floating);
                SameArithmetic((width * floating, height * floating), floating * actual);
                SameArithmetic((width / floating, height / floating), actual / floating);
            }
            Same((width, height), actual);
        }
    }

    [Theory]
    [InlineData(731)]
    [InlineData(42)]
    public void FloatArithmeticAndRoundingUseScalarComponents(int seed)
    {
        var random = new Random(seed);
        var sizes = new List<(float Width, float Height)>();
        foreach (float width in Floats)
        foreach (float height in Floats)
            sizes.Add((width, height));
        for (int i = 0; i < 256; i++)
        {
            float Next() => BitConverter.Int32BitsToSingle(
                (int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1));
            sizes.Add((Next(), Next()));
        }
        foreach (var expected in sizes)
        {
            var actual = new SizeFloat(expected.Width, expected.Height);
            foreach (var second in sizes)
            {
                var secondActual = new SizeFloat(second.Width, second.Height);
                var sum = (expected.Width + second.Width, expected.Height + second.Height);
                var difference = (expected.Width - second.Width, expected.Height - second.Height);
                var product = (expected.Width * second.Width, expected.Height * second.Width);
                SameArithmetic(sum, SizeFloat.Add(actual, secondActual));
                SameArithmetic(difference, SizeFloat.Subtract(actual, secondActual));
                SameArithmetic(sum, actual + secondActual);
                SameArithmetic(difference, actual - secondActual);
                SameArithmetic(product, actual * second.Width);
                SameArithmetic(product, second.Width * actual);
                SameArithmetic((expected.Width / second.Width, expected.Height / second.Width), actual / second.Width);
                bool equal = expected.Width == second.Width && expected.Height == second.Height;
                Assert.Equal(equal, actual.Equals(secondActual));
                Assert.Equal(equal, actual == secondActual);
                Assert.Equal(!equal, actual != secondActual);
            }
            Same((unchecked((int)Math.Ceiling(expected.Width)), unchecked((int)Math.Ceiling(expected.Height))), SizeInt.Ceiling(actual));
            Same((unchecked((int)Math.Round(expected.Width)), unchecked((int)Math.Round(expected.Height))), SizeInt.Round(actual));
            Same((unchecked((int)expected.Width), unchecked((int)expected.Height)), SizeInt.Truncate(actual));
            Same(expected, actual);
        }
    }

    [Theory]
    [InlineData("en-US", "{Width=-3.25, Height=7.5}")]
    [InlineData("ru-RU", "{Width=-3,25, Height=7,5}")]
    public void FormattingUsesCurrentCulture(string culture, string expected)
    {
        var previous = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo(culture);
            Assert.Equal("{Width=-3, Height=7}", new SizeInt(-3, 7).ToString());
            Assert.Equal(expected, new SizeFloat(-3.25f, 7.5f).ToString());
        }
        finally { CultureInfo.CurrentCulture = previous; }
    }
}
