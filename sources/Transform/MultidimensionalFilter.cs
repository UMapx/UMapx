using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the multidimensional filter.
    /// </summary>
    [Serializable]
    public class MultidimensionalFilter
    {
        #region Initialize
        /// <summary>
        /// Initializes the multidimensional filter.
        /// </summary>
        /// <param name="filter">IFilter</param>
        public MultidimensionalFilter(IFilter filter)
        {
            this.Filter = filter;
        }
        /// <summary>
        /// Gets or sets filter.
        /// </summary>
        public IFilter Filter { get; set; }
        #endregion

        #region Transform methods
        /// <summary>
        /// Forward multidimensional filter.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <returns>Jagged array</returns>
        public void Apply(params float[][] A)
        {
            int count = A.Length;

            for (int i = 0; i < count; i++)
            {
                Filter.Apply(A[i]);
            }
        }
        /// <summary>
        /// Forward multidimensional filter.
        /// </summary>
        /// <param name="A">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public void Apply(params float[][,] A)
        {
            int count = A.Length;

            for (int i = 0; i < count; i++)
            {
                Filter.Apply(A[i]);
            }
        }
        /// <summary>
        /// Forward multidimensional filter.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <returns>Jagged array</returns>
        public void Apply(params ComplexF[][] A)
        {
            int count = A.Length;

            for (int i = 0; i < count; i++)
            {
                Filter.Apply(A[i]);
            }
        }
        /// <summary>
        /// Forward multidimensional filter.
        /// </summary>
        /// <param name="A">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public void Apply(params ComplexF[][,] A)
        {
            int count = A.Length;

            for (int i = 0; i < count; i++)
            {
                Filter.Apply(A[i]);
            }
        }
        #endregion
    }
}
