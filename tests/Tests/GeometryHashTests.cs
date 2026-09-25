using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class GeometryHashTests
{
    private static void EqualKeys<T>(T first, T second) where T : struct, IEquatable<T>
    {
        Assert.True(first.Equals(second));
        Assert.Equal(first.GetHashCode(), second.GetHashCode());
        var dictionary = new Dictionary<T, string> { [first] = "value" };
        Assert.Equal("value", dictionary[second]);
        dictionary[second] = "replacement";
        Assert.Single(dictionary);
        Assert.Equal("replacement", dictionary[first]);
    }

    [Fact]
    public void EqualIntegerGeometryValuesWorkAsDictionaryKeysIncludingExtremes()
    {
        int[] values = { int.MinValue, int.MaxValue, -65536, -1, 0, 1, 65536 };
        foreach (int x in values)
        foreach (int y in values)
        {
            EqualKeys(new PointInt(x, y), new PointInt(x, y));
            EqualKeys(new SizeInt(x, y), new SizeInt(x, y));
            EqualKeys(new RectangleInt(x, y, y, x), new RectangleInt(x, y, y, x));
        }
    }

    [Fact]
    public void EqualFloatGeometryValuesWorkAsDictionaryKeysIncludingExtremes()
    {
        float[] values = { float.MinValue, float.MaxValue, -float.Epsilon, float.Epsilon,
            float.NegativeInfinity, float.PositiveInfinity, -3.5f, 2.25f, 0 };
        foreach (float x in values)
        foreach (float y in values)
        {
            EqualKeys(new PointFloat(x, y), new PointFloat(x, y));
            EqualKeys(new SizeFloat(x, y), new SizeFloat(x, y));
            EqualKeys(new RectangleFloat(x, y, y, x), new RectangleFloat(x, y, y, x));
        }
    }

    [Fact]
    public void SignedZerosAreInterchangeableInEveryComponentOfDictionaryKeys()
    {
        float negativeZero = BitConverter.Int32BitsToSingle(int.MinValue);
        EqualKeys(new PointFloat(0, 3.5f), new PointFloat(negativeZero, 3.5f));
        EqualKeys(new PointFloat(3.5f, 0), new PointFloat(3.5f, negativeZero));
        EqualKeys(new SizeFloat(0, 3.5f), new SizeFloat(negativeZero, 3.5f));
        EqualKeys(new SizeFloat(3.5f, 0), new SizeFloat(3.5f, negativeZero));
        EqualKeys(new RectangleFloat(0, 2, 3, 4), new RectangleFloat(negativeZero, 2, 3, 4));
        EqualKeys(new RectangleFloat(1, 0, 3, 4), new RectangleFloat(1, negativeZero, 3, 4));
        EqualKeys(new RectangleFloat(1, 2, 0, 4), new RectangleFloat(1, 2, negativeZero, 4));
        EqualKeys(new RectangleFloat(1, 2, 3, 0), new RectangleFloat(1, 2, 3, negativeZero));
    }

    [Theory]
    [InlineData(typeof(PointInt), 2, false)]
    [InlineData(typeof(PointFloat), 2, true)]
    [InlineData(typeof(SizeInt), 2, false)]
    [InlineData(typeof(SizeFloat), 2, true)]
    [InlineData(typeof(RectangleInt), 4, false)]
    [InlineData(typeof(RectangleFloat), 4, true)]
    public void HashDistinguishesUnitChangesInEachComponentAndTheirOrder(Type type, int count, bool floating)
    {
        object Zero() => floating ? (object)0f : 0;
        object One() => floating ? (object)1f : 1;
        var components = Enumerable.Range(0, count).Select(_ => Zero()).ToArray();
        int originalHash = Activator.CreateInstance(type, components)!.GetHashCode();
        var hashes = new HashSet<int> { originalHash };
        for (int i = 0; i < count; i++)
        {
            components[i] = One();
            int changedHash = Activator.CreateInstance(type, components)!.GetHashCode();
            Assert.True(hashes.Add(changedHash), $"{type.Name}, component {i}");
            components[i] = Zero();
        }
    }

    [Fact]
    public void HashesOfNaNComponentsAreRepeatableWithoutChangingEquality()
    {
        float[] values = { float.NaN, BitConverter.Int32BitsToSingle(0x7fc00001),
            BitConverter.Int32BitsToSingle(unchecked((int)0xffc12345)) };
        foreach (float value in values)
        {
            var point = new PointFloat(value, 3);
            Assert.Equal(point.GetHashCode(), point.Clone().GetHashCode());
            Assert.False(point.Equals(point.Clone()));
            var size = new SizeFloat(3, value);
            Assert.Equal(size.GetHashCode(), size.Clone().GetHashCode());
            Assert.False(size.Equals(size.Clone()));
            var rectangle = new RectangleFloat(1, 2, value, 4);
            Assert.Equal(rectangle.GetHashCode(), rectangle.Clone().GetHashCode());
            Assert.False(rectangle.Equals(rectangle.Clone()));
        }
    }
}
