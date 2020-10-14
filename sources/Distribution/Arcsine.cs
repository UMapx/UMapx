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
        public double Mean
        {
            get { return 0.5; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { return 0.5; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return 1.0 / 8; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return double.NaN; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get { return -1.5; }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(0, 1); }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        public double Entropy
        {
            get { return Math.Log(Math.PI / 4); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            return 2.0 / Math.PI * Math.Asin(Math.Sqrt(x));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 1 && x > 0)
            {
                return 1.0 / (Math.PI * Math.Sqrt(x * (1 - x)));
            }
            return double.NaN;
        }
        #endregion
    }
}
