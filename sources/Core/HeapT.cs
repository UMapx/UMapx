using System;
using System.Collections.Generic;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a heap.
    /// </summary>
    /// <typeparam name="T">Type</typeparam>
    [Serializable]
    public class Heap<T>
    {
        #region Private data
        private readonly List<T> data = new List<T>();
        private readonly IComparer<T> comparer;
        #endregion

        #region Constructor
        /// <summary>
        /// Initializes the heap.
        /// </summary>
        /// <param name="comparer">Comparer</param>
        public Heap(IComparer<T> comparer) => this.comparer = comparer;
        #endregion

        #region Properties and methods
        /// <summary>
        /// Gets the heap items count.
        /// </summary>
        public int Count => data.Count;

        /// <summary>
        /// Adds an item to the heap.
        /// </summary>
        /// <param name="item">Item</param>
        public void Add(T item)
        {
            data.Add(item);
            int ci = data.Count - 1;

            while (ci > 0)
            {
                int pi = (ci - 1) / 2;
                if (comparer.Compare(data[ci], data[pi]) >= 0) break;
                (data[ci], data[pi]) = (data[pi], data[ci]);
                ci = pi;
            }
        }

        /// <summary>
        /// Returns the element at the top of the heap without removing it.
        /// </summary>
        /// <returns>Element at the top of the heap</returns>
        public T Peek()
        {
            if (Count == 0)
                throw new InvalidOperationException("Cannot peek from an empty heap");
            return data[0];
        }

        /// <summary>
        /// Extracts an item from the heap.
        /// </summary>
        /// <returns>The removed root element</returns>
        public T Extract()
        {
            if (Count == 0)
                throw new InvalidOperationException("Cannot extract from an empty heap");

            int li = data.Count - 1;
            T front = data[0];
            data[0] = data[li];
            data.RemoveAt(li);
            li--;
            int pi = 0;

            while (true)
            {
                int ci = pi * 2 + 1;
                if (ci > li) break;
                int rc = ci + 1;
                if (rc <= li && comparer.Compare(data[rc], data[ci]) < 0)
                    ci = rc;
                if (comparer.Compare(data[pi], data[ci]) <= 0) break;
                (data[pi], data[ci]) = (data[ci], data[pi]);
                pi = ci;
            }

            return front;
        }
        #endregion
    }
}
