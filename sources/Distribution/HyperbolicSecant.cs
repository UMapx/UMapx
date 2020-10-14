using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the hyperbolic secant distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    /// </remarks>
    /// </summary>
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
        public double Mean
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return 1; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return 0; }
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
            get { return 2; }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(Double.NegativeInfinity, Double.PositiveInfinity); }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        public double Entropy
        {
            get { return (4.0 / Math.PI) * Maths.G; }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            double angle = Math.Atan(Math.Exp(x * Math.PI / 2.0));
            return 2 * angle / Math.PI;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return 0.5 * Maths.Sch(x * (Math.PI / 2.0));
        }
        #endregion
    }
}
