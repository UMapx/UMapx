using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the uniform distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Uniform : IDistribution
    {
        #region Private data
        private double a = 1;
        private double b = 1;
        #endregion

        #region Uniform components
        /// <summary>
        /// Initializes the uniform distribution.
        /// </summary>
        public Uniform() { }
        /// <summary>
        /// Initializes the uniform distribution.
        /// </summary>
        /// <param name="a">Shift parameter a</param>
        /// <param name="b">Shift parameter b</param>
        public Uniform(double a, double b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the shift parameter a.
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
        /// Gets or sets the shift parameter b.
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
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return (a + b) / 2.0;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return Math.Pow(b - a, 2) / 12.0;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return (a - b) / 2.0;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Mean;
            }
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
                return -1.2;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < a)
            {
                return 0;
            }
            else if (x > b)
            {
                return 0;
            }
            return 1.0 / (b - a);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < a)
            {
                return 0;
            }
            else if (x > b)
            {
                return 1;
            }
            return (x - a) / (b - a);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Math.Log(b - a);
            }
        }
        #endregion
    }
}
