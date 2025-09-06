using System;
using System.Collections.Generic;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a heap set.
    /// </summary>
    /// <typeparam name="T">Type</typeparam>
    [Serializable]
    public class HeapSet<T>
    {
        #region Private data
        private readonly Heap<T> low;  // max-heap
        private readonly Heap<T> high; // min-heap
        private readonly IComparer<T> comparer;
        #endregion

        #region Constructor
        /// <summary>
        /// Initializes the heap set.
        /// </summary>
        /// <param name="comparer">Comparer</param>
        public HeapSet(IComparer<T> comparer)
        {
            this.comparer = comparer;
            low = new Heap<T>(Comparer<T>.Create((a, b) => comparer.Compare(b, a)));
            high = new Heap<T>(comparer);
        }
        #endregion

        #region Properties and methods
        /// <summary>
        /// Gets the heap set items count.
        /// </summary>
        public int Count => low.Count + high.Count;

        /// <summary>
        /// Adds an item to the heap set.
        /// </summary>
        /// <param name="item">Item</param>
        public void Add(T item)
        {
            if (low.Count == 0 || comparer.Compare(item, low.Peek()) <= 0)
                low.Add(item);
            else
                high.Add(item);
        }

        /// <summary>
        /// Balances the heap set.
        /// </summary>
        /// <param name="targetLowSize">Target low size</param>
        public void Balance(int targetLowSize)
        {
            while (low.Count > targetLowSize)
                high.Add(low.Extract());

            while (low.Count < targetLowSize && high.Count > 0)
                low.Add(high.Extract());
        }

        /// <summary>
        /// Returns the rank.
        /// </summary>
        /// <returns>Element</returns>
        public T GetRank() => low.Peek();
        #endregion
    }
}
