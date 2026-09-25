using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class HeapTests
{
    [Theory]
    [InlineData(1)]
    [InlineData(31)]
    [InlineData(314159)]
    public void HeapsAndRankContainersMatchSortedMultisets(int seed)
    {
        var random = new Random(seed);
        var values = Enumerable.Range(0, 101).Select(_ => random.Next(-20, 21)).ToArray();
        var sorted = values.OrderBy(v => v).ToArray();
        foreach (bool reverse in new[]
        {
            false,
            true
        }

        )
        {
            var comparer = Comparer<int>.Create((a, b) => reverse ? b.CompareTo(a) : a.CompareTo(b));
            var heap = new Heap<int>(comparer);
            foreach (int value in values)
                heap.Add(value);
            Assert.Equal(values.Length, heap.Count);
            foreach (int value in reverse ? sorted.Reverse() : sorted)
            {
                Assert.Equal(value, heap.Peek());
                Assert.Equal(value, heap.Extract());
            }

            Assert.Equal(0, heap.Count);
            Assert.Throws<InvalidOperationException>(() => heap.Extract());
            Assert.Throws<InvalidOperationException>(() => heap.Peek());
        }

        var set = new HeapSet<int>(Comparer<int>.Default);
        foreach (int value in values)
            set.Add(value);
        Assert.Equal(values.Length, set.Count);
        foreach (int rank in new[]
        {
            0,
            50,
            100,
            1,
            99,
            30,
            70
        }

        )
        {
            set.Balance(rank + 1);
            Assert.Equal(sorted[rank], set.GetRank());
        }
    }
}
