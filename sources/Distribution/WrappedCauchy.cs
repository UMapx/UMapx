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
        private float mu;
        private float gamma;
        #endregion

        #region Wrapped distribution
        /// <summary>
        /// Initializes the compact Cauchy distribution.
        /// </summary>
        /// <param name="mu">Parameter μ</param>
        /// <param name="gamma">Parameter γ > 0</param>
        public WrappedCauchy(float mu, float gamma)
        {
            this.mu = mu;
            this.gamma = gamma;
        }
        /// <summary>
        /// Gets or sets the value of the parameter μ.
        /// </summary>
        public float Mu
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
        public float Gamma
        {
            get
            {
                return this.gamma;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.gamma = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return mu; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return 1 - (float)Math.Exp(-gamma); }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(-Maths.Pi, Maths.Pi); }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get { return (float)Math.Log(2 * Math.PI * (1 - Math.Exp(-2 * gamma))); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            float constant = (float)(1.0 / (2 * Math.PI));
            return constant * (float)Math.Sinh(gamma) / (float)(Math.Cosh(gamma) - Math.Cos(x - mu));
        }
        #endregion
    }
}
