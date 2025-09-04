using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the logarithmic Gaussian distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Log-normal_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LogGaussian : IDistribution
    {
        #region Private data
        private float sigma = 1;
        private float mu = 0;
        #endregion;

        #region GaussianLog components
        /// <summary>
        /// Initializes the logarithmic Gaussian distribution.
        /// </summary>
        public LogGaussian() { }
        /// <summary>
        /// Initializes the logarithmic Gaussian distribution.
        /// </summary>
        /// <param name="sigma">Standard deviation</param>
        /// <param name="mu">Mathematical expectation</param>
        public LogGaussian(float sigma, float mu)
        {
            Sigma = sigma;
            Mu = mu;
        }
        /// <summary>
        /// Gets or sets the standard deviation.
        /// </summary>
        public float Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the mathematical expectation.
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
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(0, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return MathsF.Exp(mu + sigma * sigma / 2);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return (MathsF.Exp(sigma * sigma) - 1) * MathsF.Exp(2 * mu + sigma * sigma);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return MathsF.Exp(mu);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return MathsF.Exp(mu - sigma * sigma);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return (float)(MathsF.Exp(sigma * sigma) + 2.0) * MathsF.Sqrt(MathsF.Exp(sigma * sigma) - 1.0f);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                return (float)MathsF.Exp(4 * sigma * sigma) + 2.0f * MathsF.Exp(3 * sigma * sigma) + 3.0f * MathsF.Exp(3 * sigma * sigma) - 6.0f;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < 0)
            {
                return 0;
            }
            return MathsF.Exp(MathsF.Pow((MathsF.Log(x) - mu), 2) / (-2.0f * sigma * sigma)) / (MathsF.Sqrt(2.0f * MathsF.Pi) * sigma * x);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 0.5f + 0.5f * Special.Erf((MathsF.Log(x) - mu) / MathsF.Sqrt(sigma * MathsF.Sqrt2));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return 0.5f + 0.5f * MathsF.Log(2 * MathsF.Pi * sigma * sigma) + mu;
            }
        }
        #endregion
    }
}
