using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the compact Cauchy distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wrapped_Cauchy_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class WrappedCauchy : IDistribution
    {
        #region Private data
        private double mu;
        private double gamma;
        #endregion

        #region Wrapped distribution
        /// <summary>
        /// Initializes the compact Cauchy distribution.
        /// </summary>
        /// <param name="mu">Parameter μ</param>
        /// <param name="gamma">Parameter γ > 0</param>
        public WrappedCauchy(double mu, double gamma)
        {
            this.mu = mu;
            this.gamma = gamma;
        }
        /// <summary>
        /// Gets or sets the value of the parameter μ.
        /// </summary>
        public double Mu
        {
            get
            {
                return this.mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter γ > 0.
        /// </summary>
        public double Gamma
        {
            get
            {
                return this.gamma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.gamma = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return mu; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return 1 - Math.Exp(-gamma); }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(-Math.PI, Math.PI); }
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
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { return Math.Log(2 * Math.PI * (1 - Math.Exp(-2 * gamma))); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double constant = (1.0 / (2 * Math.PI));
            return constant * Math.Sinh(gamma) / (Math.Cosh(gamma) - Math.Cos(x - mu));
        }
        #endregion
    }
}
