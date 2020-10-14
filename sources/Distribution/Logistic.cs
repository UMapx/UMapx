using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the logistic distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Logistic_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Logistic : IDistribution
    {
        #region Private data
        private double mu = 5;
        private double s = 2;
        #endregion

        #region Logistic components
        /// <summary>
        /// Initializes the logistic distribution.
        /// </summary>
        /// <param name="mu">Parameter μ</param>
        /// <param name="s">Parameter s (0, +inf]</param>
        public Logistic(double mu, double s)
        {
            Mu = mu; S = s;
        }
        /// <summary>
        /// Initializes the logistic distribution.
        /// </summary>
        public Logistic() { }
        /// <summary>
        /// Gets or sets the value of the parameter μ.
        /// </summary>
        public double Mu
        {
            get
            {
                return mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter s (0, +inf].
        /// </summary>
        public double S
        {
            get
            {
                return s;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.s = value;
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
                return mu;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return (s * s * Math.PI * Math.PI) / 3.0; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return mu;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return mu; }
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
                return 1.2;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                return Math.Log(s) + 2;
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            double z = (x - mu) / s;
            return 1.0 / (1 + Math.Exp(-z));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double z = (x - mu) / s;
            double num = Math.Exp(-z);
            double a = (1 + num);
            double den = s * a * a;

            return num / den;
        }
        #endregion
    }
}
