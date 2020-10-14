using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Gamma-distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gamma_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Gamma : IDistribution
    {
        #region Private data
        private double thetta = 1;
        private double k = 1;
        #endregion

        #region Gamma components
        /// <summary>
        /// Initializes the Gamma-distribution.
        /// </summary>
        public Gamma() { }
        /// <summary>
        /// Initializes the Gamma-distribution.
        /// </summary>
        /// <param name="thetta">Parameter θ (0, +inf)</param>
        /// <param name="k">Parameter k (0, +inf)</param>
        public Gamma(double thetta, double k)
        {
            Thetta = thetta; K = k;
        }
        /// <summary>
        /// Gets or sets the parameter θ (0, +inf).
        /// </summary>
        public double Thetta
        {
            get
            {
                return this.thetta;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.thetta = value;
            }
        }
        /// <summary>
        /// Gets or sets the parameter k (0, +inf).
        /// </summary>
        public double K
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
                return thetta * k;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return k * thetta * thetta;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (k >= 1)
                {
                    return (k - 1) * thetta;
                }
                return double.NaN;
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
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2 / Math.Sqrt(k);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6.0 / k;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Math.Pow(x, k - 1) * Math.Exp(-x / thetta) / (Special.Gamma(k) * Math.Pow(thetta, k));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Special.GammaP(k, x / thetta);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return k * thetta + (1 - k) * Math.Log(thetta) + Special.GammaLog(k); // + (1 - k) * Special.Ksi(k);
            }
        }
        #endregion
    }
}
