using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Burr distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Burr_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Burr : IDistribution
    {
        #region Private data
        private double c;
        private double k;
        #endregion

        #region Burr distribution
        /// <summary>
        /// Initializes the Burr distribution.
        /// </summary>
        /// <param name="c">Form parameter c > 0</param>
        /// <param name="k">Scale parameter k > 0</param>
        public Burr(double c, double k)
        {
            C = c; K = k;
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter c > 0.
        /// </summary>
        public double C
        {
            get
            {
                return this.c;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.c = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter k > 0.
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return k * Special.Beta(k - 1.0 / c, 1.0 + 1.0 / c); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return Math.Pow((c - 1) / (k * c + 1), 1.0 / c);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Pow(Math.Pow(2, 1.0 / k) - 1.0, 1.0 / c);
            }
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
        /// Returns the value of differential entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(double.Epsilon, Double.PositiveInfinity); }
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
                return double.NaN;
            }

            return 1.0 - Math.Pow(1.0 + Math.Pow(x, c), -k);
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
                return double.NaN;
            }

            double a = c * k;
            double b = Math.Pow(x, c - 1);
            double d = 1 + Math.Pow(x, c);
            return a * b / Math.Pow(d, k + 1);
        }
        #endregion
    }
}
