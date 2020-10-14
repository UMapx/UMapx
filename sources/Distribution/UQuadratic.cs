using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the U-quadratic distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/U-quadratic_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class UQuadratic : IDistribution
    {
        #region Private data
        double a;
        double b;
        double alpha;
        double beta;
        #endregion

        #region UQuadratic components
        /// <summary>
        /// Initializes the U-quadratic distribution.
        /// </summary>
        /// <param name="a">Parameter a ∈ (0, +inf)</param>
        /// <param name="b">Parameter b ∈ (a, +inf)</param>
        public UQuadratic(double a, double b)
        {
            A = a; B = b;
            this.alpha = 12 / Math.Pow(b - a, 3);
            this.beta = (b + a) / 2;
        }
        /// <summary>
        /// Gets or sets the value of the parameter a ∈ (0, +inf).
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter b ∈ (a, +inf).
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (b < a)
                    throw new Exception("The value of parameter b must be either greater than or equal to a");

                this.b = value;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return (int)a & (int)b;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return (a + b) / 2.0d; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { return (a + b) / 2.0d; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return (3.0d / 20.0) * Math.Pow(b - a, 2.0); }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 3.0 / 112.0 * Math.Pow(b - a, 4.0);
            }
        }
        /// <summary>
        /// Gets the value of entropy.
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
            get { return new RangeDouble(a, b); }
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

            if (x > b)
                return 1;

            return (alpha / 3) * (Math.Pow(x - beta, 3) + Math.Pow(beta - a, 3));
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

            if (x > b)
                return 0;

            return alpha * Math.Pow(x - beta, 2);
        }
        #endregion
    }
}
