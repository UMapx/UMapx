using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the triangular distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Triangular_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Triangular : IDistribution
    {
        #region Private data
        private double a;
        private double b;
        private double c;
        #endregion

        #region Triangular components
        /// <summary>
        /// Initializes the triangular distribution.
        /// </summary>
        public Triangular() { }
        /// <summary>
        /// Initializes the triangular distribution.
        /// </summary>
        /// <param name="a">Parameter a ∈ (-inf, +inf)</param>
        /// <param name="b">Parameter b ∈ (-inf, +inf)</param>
        /// <param name="c">Parameter c ∈ (-inf, +inf)</param>
        public Triangular(double a, double b, double c)
        {
            A = a; B = b; C = c;
        }
        /// <summary>
        /// Gets or sets the value of the parameter a ∈ (-inf, +inf).
        /// </summary>
        public double A
        {
            get
            {
                return a;
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter b ∈ (-inf, +inf).
        /// </summary>
        public double B
        {
            get
            {
                return b;
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter c ∈ (-inf, +inf).
        /// </summary>
        public double C
        {
            get
            {
                return c;
            }
            set
            {
                this.c = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(a, b);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return (a + b + c) / 3.0; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return (a * a + b * b + c * c - a * b - a * c - b * c) / 18; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                double median;
                if (c >= (a + b) / 2.0)
                {
                    median = a + Math.Sqrt((b - a) * (c - a)) / 1.4142135623730950488016887242097;
                }
                else
                {
                    median = b - Math.Sqrt((b - a) * (b - c)) / 1.4142135623730950488016887242097;
                }

                return median;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return c; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                double k1 = (a + b - 2 * c) * (2 * a - b - c) * (a - 2 * b + c);
                double k2 = 5 * (a * a + b * b + c * c - a * b - a * c - b * c);
                return 1.4142135623730950488016887242097 * k1 / Math.Pow(k2, 1.5);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return -0.6;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                return 0.5 + Math.Log((b - a) / 2);
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < a)
                return 0;

            if (x >= a && x <= c)
                return ((x - a) * (x - a)) / ((b - a) * (c - a));

            if (x > c && x <= b)
                return 1 - ((b - x) * (b - x)) / ((b - a) * (b - c));

            return 1;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < a)
                return 0;

            if (x >= a && x <= c)
                return (2 * (x - a)) / ((b - a) * (c - a));

            if (x > c && x <= b)
                return (2 * (b - x)) / ((b - a) * (b - c));

            return 0;
        }
        #endregion
    }
}
