using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the hyperbolic secant distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    /// </remarks>
    [Serializable]
    public class HyperbolicSecant : IDistribution
    {
        #region Distribution
        /// <summary>
        /// Initializes the hyperbolic secant distribution.
        /// </summary>
        public HyperbolicSecant() { }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        /// <remarks>
        /// For the standard hyperbolic secant distribution, the variance equals 1.
        /// </remarks>
        public float Variance
        {
            get { return 1f; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        public float Excess
        {
            get { return -1f; }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity); }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <remarks>
        /// This value corresponds to the hyperbolic secant distribution with scale σ = 1.
        /// </remarks>
        public float Entropy
        {
            get { return Maths.Log(4f); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            float angle = Maths.Atan(Maths.Exp(x * Maths.Pi / 2.0f));
            return 2 * angle / Maths.Pi;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            return 0.5f * Maths.Sech(x * (Maths.Pi / 2.0f));
        }
        #endregion
    }
}
