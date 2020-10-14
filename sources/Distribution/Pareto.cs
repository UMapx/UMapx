using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Pareto distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Pareto_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Pareto : IDistribution
    {
        #region Private data
        private double xm = 1;
        private double k = 1;
        #endregion

        #region Pareto components
        /// <summary>
        /// Initializes the Pareto distribution.
        /// </summary>
        public Pareto() { }
        /// <summary>
        /// Initializes the Pareto distribution.
        /// </summary>
        /// <param name="xm">Scale factor θ (0, +inf)</param>
        /// <param name="k">Parameter k (0, +inf)</param>
        public Pareto(double xm, double k)
        {
            Xm = xm; K = k;
        }
        /// <summary>
        /// Gets or sets the scale factor Xm (0, +inf).
        /// </summary>
        public double Xm
        {
            get
            {
                return this.xm;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.xm = value;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor k (0, +inf).
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
                return (k * xm) / (k - 1);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                if (k > 2)
                {
                    return Maths.Pow(xm / k - 1) * (k / (k - 2));
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return xm;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return xm * Maths.Sqrt(2, k);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                if (k > 3)
                {
                    return 2 * (1 + k) / (k - 3) * Math.Sqrt((k - 2) / k);
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                double k2 = k * k;
                double k3 = k2 * k;
                return 6 * (k3 + k2 + 6 * k - 2) / (k * (k - 3) * (k - 4));
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < xm)
            {
                return 0;
            }
            return k * Math.Pow(xm, k) / Math.Pow(x, k + 1);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < xm)
            {
                return 0;
            }
            return 1.0 - Math.Pow(xm / x, k);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return k * xm + (1 - k) * Math.Log(xm) + Special.GammaLog(k); // + (1 - k) * Special.Ksi(k);
            }
        }
        #endregion
    }
}
