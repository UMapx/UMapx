﻿using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the filter interface.
    /// </summary>
    public interface IFilter
    {
        #region Interface
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        void Apply(float[] data);
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        void Apply(float[,] data);
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        void Apply(Complex32[] data);
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        void Apply(Complex32[,] data);
        #endregion
    }
}
