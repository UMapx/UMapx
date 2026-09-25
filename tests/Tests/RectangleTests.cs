using System.Globalization;
using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class RectangleTests
{
    [Fact]
    public void ClonesPreserveAllComponentsAndCanBeModifiedIndependently()
    {
        foreach (var expected in IntegerCases())
        {
            var original = Convert(expected);
            var clone = original.Clone();
            Same(expected, clone);
            Same(expected, (RectangleInt)((ICloneable)original).Clone());
            clone.X = 123;
            clone.Height = -456;
            Same(expected, original);
        }
        foreach (var expected in FloatCases())
        {
            var original = Convert(expected);
            var clone = original.Clone();
            Same(expected, clone);
            Same(expected, (RectangleFloat)((ICloneable)original).Clone());
            clone.X = 123;
            clone.Height = -456;
            Same(expected, original);
        }
    }

    private static void Same((int X, int Y, int Width, int Height) expected, RectangleInt actual)
    {
        Assert.Equal(expected.X, actual.X);
        Assert.Equal(expected.Y, actual.Y);
        Assert.Equal(expected.Width, actual.Width);
        Assert.Equal(expected.Height, actual.Height);
    }

    private static void Same((float X, float Y, float Width, float Height) expected, RectangleFloat actual)
    {
        SameBits(expected.X, actual.X);
        SameBits(expected.Y, actual.Y);
        SameBits(expected.Width, actual.Width);
        SameBits(expected.Height, actual.Height);
    }

    private static (int X, int Y, int Width, int Height)[] IntegerCases(int seed = 731)
    {
        var cases = new List<(int X, int Y, int Width, int Height)>
        {
            (0, 0, 0, 0), (3, 4, 0, 0), (0, 0, 10, 20), (10, 20, 5, 5),
            (10, 0, 4, 3), (-4, -3, 10, 20), (2, 3, -4, 5), (2, 3, 4, -5),
            (2, 3, -4, -5), (5, 5, 0, 4), (5, 5, 4, 0),
            (int.MaxValue, int.MinValue, 4, -5), (int.MinValue, int.MaxValue, -4, 5),
            (-1, -2, int.MaxValue, int.MinValue)
        };
        var random = new Random(seed);
        for (int i = 0; i < 60; i++)
            cases.Add((random.Next(-50, 50), random.Next(-50, 50),
                random.Next(-10, 70), random.Next(-10, 70)));
        for (int i = 0; i < 64; i++)
            cases.Add(((int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1),
                (int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1),
                (int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1),
                (int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1)));
        return cases.ToArray();
    }

    private static (float X, float Y, float Width, float Height)[] FloatCases(int seed = 731)
    {
        var cases = IntegerCases(seed).Select(r => (X: (float)r.X, Y: (float)r.Y, Width: (float)r.Width, Height: (float)r.Height)).ToList();
        cases.AddRange(new (float X, float Y, float Width, float Height)[]
        {
            (-2.75f, 3.125f, 4.5f, 1.75f),
            (-0.0f, 0, 0, -0.0f),
            (float.Epsilon, -float.Epsilon, float.Epsilon, float.Epsilon),
            (float.MaxValue, float.MinValue, float.MaxValue, float.MaxValue),
            (float.PositiveInfinity, 0, 2, 3),
            (0, float.NegativeInfinity, 2, 3),
            (0, 0, float.PositiveInfinity, 3),
            (float.NaN, 0, 2, 3),
            (0, float.NaN, 2, 3),
            (0, 0, float.NaN, 3),
            (0, 0, 2, float.NaN)
        });
        var random = new Random(seed);
        for (int i = 0; i < 64; i++)
        {
            float Next() => BitConverter.Int32BitsToSingle((int)random.NextInt64(int.MinValue, (long)int.MaxValue + 1));
            cases.Add((Next(), Next(), Next(), Next()));
        }
        return cases.ToArray();
    }

    [Fact]
    public void IntegerConstructionPropertiesAndEqualityPreserveComponents()
    {
        foreach (var expected in IntegerCases())
        {
            var actual = new RectangleInt(new PointInt(expected.X, expected.Y),
                new SizeInt(expected.Width, expected.Height));
            Same(expected, actual);
            Assert.Equal(expected.X, actual.Left);
            Assert.Equal(expected.Y, actual.Top);
            Assert.Equal(unchecked(expected.X + expected.Width), actual.Right);
            Assert.Equal(unchecked(expected.Y + expected.Height), actual.Bottom);
            Assert.Equal(new PointInt(expected.X, expected.Y), actual.Location);
            Assert.Equal(new SizeInt(expected.Width, expected.Height), actual.Size);
            Same((expected.X, expected.Y, unchecked(expected.Width - expected.X), unchecked(expected.Height - expected.Y)),
                RectangleInt.FromLTRB(expected.X, expected.Y, expected.Width, expected.Height));
            Assert.True(actual.Equals((object)actual));
            Assert.True(((IEquatable<RectangleInt>)actual).Equals(actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals((object)expected));

            var changed = actual;
            changed.Location = new PointInt(7, -11);
            changed.Size = new SizeInt(13, 17);
            Same((7, -11, 13, 17), changed);
            changed.X = -9; changed.Y = 3; changed.Width = 0; changed.Height = -4;
            Same((-9, 3, 0, -4), changed);
            Same(expected, actual);
        }
    }

    [Fact]
    public void FloatConstructionPropertiesAndEqualityPreserveComponents()
    {
        foreach (var expected in FloatCases())
        {
            var actual = new RectangleFloat(new PointFloat(expected.X, expected.Y),
                new SizeFloat(expected.Width, expected.Height));
            Same(expected, actual);
            Assert.Equal(expected.X, actual.Left);
            Assert.Equal(expected.Y, actual.Top);
            Assert.Equal(unchecked(expected.X + expected.Width), actual.Right);
            Assert.Equal(unchecked(expected.Y + expected.Height), actual.Bottom);
            Assert.Equal(expected.Width <= 0 || expected.Height <= 0, actual.IsEmpty);
            Assert.Equal(expected.X, actual.Location.X);
            Assert.Equal(expected.Y, actual.Location.Y);
            Assert.Equal(expected.Width, actual.Size.Width);
            Assert.Equal(expected.Height, actual.Size.Height);
            Same((expected.X, expected.Y, expected.Width - expected.X, expected.Height - expected.Y),
                RectangleFloat.FromLTRB(expected.X, expected.Y, expected.Width, expected.Height));
            bool reflexive = !float.IsNaN(expected.X) && !float.IsNaN(expected.Y)
                && !float.IsNaN(expected.Width) && !float.IsNaN(expected.Height);
            Assert.Equal(reflexive, actual.Equals((object)actual));
            Assert.Equal(reflexive, ((IEquatable<RectangleFloat>)actual).Equals(actual));
            Assert.False(actual.Equals(null));
            Assert.False(actual.Equals((object)expected));

            var changed = actual;
            changed.Location = new PointFloat(7.25f, -11.5f);
            changed.Size = new SizeFloat(13.75f, 17.5f);
            Same((7.25f, -11.5f, 13.75f, 17.5f), changed);
            changed.X = -9.5f; changed.Y = 3.25f; changed.Width = 0; changed.Height = -4;
            Same((-9.5f, 3.25f, 0, -4), changed);
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
    public void IntegerOffsetAndInflationWrapOverflowInEachComponent()
    {
        foreach (var rectangle in IntegerCases())
        foreach (int amount in new[] { -20, 0, 7, int.MinValue, int.MaxValue })
        {
            RectangleInt original = Convert(rectangle);
            var expected = (X: unchecked(rectangle.X - amount), Y: unchecked(rectangle.Y + 3),
                Width: unchecked(rectangle.Width + 2 * amount), Height: unchecked(rectangle.Height - 6));
            Same(expected, RectangleInt.Inflate(original, amount, -3));
            var actual = original;
            actual.Inflate(new SizeInt(amount, -3));
            Same(expected, actual);
            expected = (unchecked(expected.X + amount), unchecked(expected.Y + 2), expected.Width, expected.Height);
            actual.Offset(new PointInt(amount, 2));
            Same(expected, actual);
            expected = (unchecked(expected.X + 4), unchecked(expected.Y + amount), expected.Width, expected.Height);
            actual.Offset(4, amount);
            Same(expected, actual);
            Same(rectangle, original);
        }
    }

    [Fact]
    public void FloatOffsetAndInflationUseEachComponent()
    {
        foreach (var rectangle in FloatCases())
        foreach (float amount in new[] { -20.25f, 0, 7.5f, float.MaxValue, float.NaN })
        {
            RectangleFloat original = Convert(rectangle);
            var expected = (X: rectangle.X - amount, Y: rectangle.Y + 3.5f,
                Width: rectangle.Width + 2 * amount, Height: rectangle.Height - 7);
            SameArithmetic(expected, RectangleFloat.Inflate(original, amount, -3.5f));
            var actual = original;
            actual.Inflate(new SizeFloat(amount, -3.5f));
            SameArithmetic(expected, actual);
            expected = (expected.X + amount, expected.Y + 2.25f, expected.Width, expected.Height);
            actual.Offset(new PointFloat(amount, 2.25f));
            SameArithmetic(expected, actual);
            expected = (expected.X + 4.5f, expected.Y + amount, expected.Width, expected.Height);
            actual.Offset(4.5f, amount);
            SameArithmetic(expected, actual);
            Same(rectangle, original);
        }
    }

    [Theory]
    [InlineData(-2.5f, 3.5f, 4.5f, -5.5f)]
    [InlineData(-2.9f, 3.1f, 0.1f, -0.1f)]
    [InlineData(0, 0, 0, 0)]
    [InlineData(16777215, -16777215, 100000, 200000)]
    public void FloatToIntegerRoundingUsesScalarConversions(float x, float y, float width, float height)
    {
        var actual = new RectangleFloat(x, y, width, height);
        Same((unchecked((int)Math.Ceiling(x)), unchecked((int)Math.Ceiling(y)), unchecked((int)Math.Ceiling(width)), unchecked((int)Math.Ceiling(height))), RectangleInt.Ceiling(actual));
        Same((unchecked((int)Math.Round(x)), unchecked((int)Math.Round(y)), unchecked((int)Math.Round(width)), unchecked((int)Math.Round(height))), RectangleInt.Round(actual));
        Same((unchecked((int)x), unchecked((int)y), unchecked((int)width), unchecked((int)height)), RectangleInt.Truncate(actual));
    }

    [Fact]
    public void ConversionsPreserveComponents()
    {
        foreach (var rectangle in IntegerCases())
        {
            RectangleInt integer = Convert(rectangle);
            RectangleFloat floating = integer;
            Same(((float)rectangle.X, (float)rectangle.Y, (float)rectangle.Width, (float)rectangle.Height), floating);
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
    [InlineData("en-US", "{X=-3.25,Y=4.5,Width=5.75,Height=6.125}")]
    [InlineData("ru-RU", "{X=-3,25,Y=4,5,Width=5,75,Height=6,125}")]
    public void StringRepresentationUsesCurrentCulture(string culture, string expected)
    {
        var previous = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo(culture);
            Assert.Equal("{X=-3,Y=4,Width=5,Height=6}", new RectangleInt(-3, 4, 5, 6).ToString());
            Assert.Equal(expected,
                new RectangleFloat(-3.25f, 4.5f, 5.75f, 6.125f).ToString());
        }
        finally { CultureInfo.CurrentCulture = previous; }
    }

    private static void SameArithmetic((float X, float Y, float Width, float Height) expected, RectangleFloat actual)
    {
        static void Compare(float expected, float actual)
        {
            if (float.IsNaN(expected)) Assert.True(float.IsNaN(actual));
            else SameBits(expected, actual);
        }
        Compare(expected.X, actual.X);
        Compare(expected.Y, actual.Y);
        Compare(expected.Width, actual.Width);
        Compare(expected.Height, actual.Height);
    }

    private static RectangleInt Convert((int X, int Y, int Width, int Height) r) => new RectangleInt(r.X, r.Y, r.Width, r.Height);
    private static RectangleFloat Convert((float X, float Y, float Width, float Height) r) => new RectangleFloat(r.X, r.Y, r.Width, r.Height);

    private static void SameBits(float expected, float actual) =>
        Assert.Equal(BitConverter.SingleToInt32Bits(expected), BitConverter.SingleToInt32Bits(actual));

    [Fact]
    public void RoundingExtremeFloatComponentsUsesRuntimeNumericConversions()
    {
        foreach (var rectangle in FloatCases())
        {
            var actual = Convert(rectangle);
            Same((unchecked((int)Math.Ceiling(rectangle.X)), unchecked((int)Math.Ceiling(rectangle.Y)), unchecked((int)Math.Ceiling(rectangle.Width)), unchecked((int)Math.Ceiling(rectangle.Height))), RectangleInt.Ceiling(actual));
            Same((unchecked((int)Math.Round(rectangle.X)), unchecked((int)Math.Round(rectangle.Y)), unchecked((int)Math.Round(rectangle.Width)), unchecked((int)Math.Round(rectangle.Height))), RectangleInt.Round(actual));
            Same((unchecked((int)rectangle.X), unchecked((int)rectangle.Y), unchecked((int)rectangle.Width), unchecked((int)rectangle.Height)), RectangleInt.Truncate(actual));
        }
    }
}
