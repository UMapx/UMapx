using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the log-logistic distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Log-logistic_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LogLogistic : IDistribution
    {
        #region Private data
        private double a = 1;
        private double b = 1;
        #endregion

        #region LogLogistic components
        /// <summary>
        /// Initializes the log-logistic distribution.
        /// </summary>
        public LogLogistic() { }
        /// <summary>
        /// Initializes the log-logistic distribution.
        /// </summary>
        /// <param name="a">Parameter a</param>
        /// <param name="b">Parameter b</param>
        public LogLogistic(double a, double b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the value of parameter a.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of parameter b.
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

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
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
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
                if (b > 1)
                {
                    return a * Math.Pow((b - 1) / (b + 1), 1 / b);
                }
                return 0;
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
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return (b / a) * Math.Pow(x / a, b - 1) / (1.0 + Math.Pow(Math.Pow(x / a, b), 2));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            return 1.0 / (1 + Math.Pow(x / a, -b));
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
