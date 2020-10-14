using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the logarithmic distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Logarithmic_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Logarithmic : IDistribution
    {
        #region Private data
        private double p = 0.66;
        #endregion

        #region Logarithmic components
        /// <summary>
        /// Initializes the logarithmic distribution.
        /// </summary>
        /// <param name="p">Parameter</param>
        public Logarithmic(double p)
        {
            P = p;
        }
        /// <summary>
        /// Gets or sets the value of the parameter p ∈ (0, 1].
        /// </summary>
        public double P
        {
            get
            {
                return p;
            }
            set
            {
                if (value <= 0 || value > 1)
                    throw new Exception("Invalid argument value");

                this.p = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(1, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return (-1 / Math.Log(1 - p)) * p / (1 - p);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                double k1 = p + Math.Log(1 - p);
                double k2 = Math.Pow(1 - p, 2) * Math.Pow(Math.Log(1 - p), 2);
                return -p * k1 / k2;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return 1; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x <= 0)
            {
                return 0;
            }
            if (x > 1)
            {
                return 0;
            }
            return 1 + Special.Beta(x + 1, 0) / Math.Log(1 - p);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x <= 0)
            {
                return 0;
            }
            if (x > 1)
            {
                return 0;
            }
            return -1 / Math.Log(1 - p) * Math.Pow(p, x) / x;
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
}
