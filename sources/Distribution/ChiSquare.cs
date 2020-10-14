using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the xi-square distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Chi-squared_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ChiSquare : IDistribution
    {
        #region Private data
        private int k = 1;
        #endregion

        #region Chi-square components
        /// <summary>
        /// Initializes the xi-square distribution.
        /// </summary>
        /// <param name="k">Degrees of freedom (0, +inf)</param>
        public ChiSquare(int k)
        {
            K = k;
        }
        /// <summary>
        /// Gets or sets the degrees of freedom (0, +inf).
        /// </summary>
        public int K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return k;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return k - 2.0 / 3;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (k >= 2)
                {
                    return k - 2;
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return 2 * k;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return Math.Sqrt(8.0 / k);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 12.0 / k;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                double k2 = k / 2.0;
                double s1 = Math.Log(2.0 * Special.Gamma(k2));
                double s2 = (1.0 - k2) * Special.DiGamma(k2);
                return k2 + s1 + s2;
            }
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
            return Special.GammaP(k / 2.0, x / 2.0);
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
                return 1;
            }
            return Special.GammaQ(k / 2.0, x / 2.0);
        }
        #endregion
    }
}
