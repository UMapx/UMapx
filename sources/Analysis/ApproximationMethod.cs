using System;

namespace UMapx.Analysis
{
    /// <summary>
    /// Approximation method.
    /// </summary>
    [Serializable]
    public enum ApproximationMethod
    {
        /// <summary>
        /// Polynomial approximation.
        /// </summary>
        Polynomial,
        /// <summary>
        /// Logarithmic approximation.
        /// </summary>
        Logarithmic,
        /// <summary>
        /// Exponential approximation.
        /// </summary>
        Exponential,
        /// <summary>
        /// Power approximation.
        /// </summary>
        Power,
    }
}
