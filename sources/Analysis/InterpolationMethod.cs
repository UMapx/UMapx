using System;

namespace UMapx.Analysis
{
    /// <summary>
    /// Interpolation method.
    /// </summary>
    [Serializable]
    public enum InterpolationMethod
    {
        /// <summary>
        /// Linear method.
        /// </summary>
        Linear,
        /// <summary>
        /// Lagrange's method.
        /// </summary>
        Lagrange,
        /// <summary>
        /// Newton's method.
        /// </summary>
        Newton,
        /// <summary>
        /// Barycentric method.
        /// </summary>
        Barycentric,
    }
}
