using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the beta distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Beta_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Beta : IDistribution
    {
        #region Private data
        private double a = 1;
        private double b = 1;
        #endregion

        #region Beta components
        /// <summary>
        /// Initializes beta distribution.
        /// </summary>
        public Beta() { }
        /// <summary>
        /// Initializes beta distribution.
        /// </summary>
        /// <param name="a">Parameter a</param>
        /// <param name="b">Parameter b</param>
        public Beta(double a, double b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the parameter a.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the parameter b.
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, 1);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return a / (a + b);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return (a * b) / Maths.Pow(a + b) / (a + b + 1);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (a > 1 && b > 1)
                {
                    return (a - 1) / (a + b - 2);
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
                return 2 * (b - a) * Math.Sqrt(a + b + 1) / (a + b + 2) / Math.Sqrt(a * b);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                double a2 = a * a, b2 = b * b, a3 = a2 * a;
                return 6 * (a3 - a2 * (2 * b - 1) + b2 * (b + 1) - 2 * a * b * (b + 2)) / (a * b * (a + b + 2) * (a + b + 3));
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x > 1)
            {
                return 0;
            }
            else if (x < 0)
            {
                return 0;
            }
            return Math.Pow(x, a - 1) * Math.Pow(1 - x, b - 1) / Special.Beta(a, b);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x > 1)
            {
                return 0;
            }
            else if (x < 0)
            {
                return 0;
            }
            return Special.Beta(a, b, x);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
}
