using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the arcsine distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Arcsine_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Arcsine : IDistribution
    {
        #region Distribution
        /// <summary>
        /// Initializes the arcsine distribution.
        /// </summary>
        public Arcsine() { }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return 0.5f; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { return 0.5f; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return 1.0f / 8; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get { return float.NaN; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get { return -1.5f; }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(0, 1); }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        public float Entropy
        {
            get { return MathF.Log(MathF.Pi / 4); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            return 2.0f / MathF.Pi * (float)Math.Asin(Math.Sqrt(x));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x > 1)
                return 1.0f;

            if (x < 0)
                return 0.0f;

            return 1.0f / (float)(Math.PI * Math.Sqrt(x * (1 - x)));
        }
        #endregion
    }
}
