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
                    throw new Exception("Invalid argument value");

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
                return Maths.Exp(mu + sigma * sigma / 2);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return (Maths.Exp(sigma * sigma) - 1) * Maths.Exp(2 * mu + sigma * sigma);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return Maths.Exp(mu);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return Maths.Exp(mu - sigma * sigma);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return (float)(Maths.Exp(sigma * sigma) + 2.0) * Maths.Sqrt(Maths.Exp(sigma * sigma) - 1.0f);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                return (float)Maths.Exp(4 * sigma * sigma) + 2.0f * Maths.Exp(3 * sigma * sigma) + 3.0f * Maths.Exp(3 * sigma * sigma) - 6.0f;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>float precision floating point number</returns>
        public float Function(float x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Maths.Exp(Maths.Pow((Maths.Log(x) - mu), 2) / (-2.0f * sigma * sigma)) / (Maths.Sqrt(2.0f * Maths.Pi) * sigma * x);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>float precision floating point number</returns>
        public float Distribution(float x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 0.5f + 0.5f * Special.Erf((Maths.Log(x) - mu) / Maths.Sqrt(sigma * Maths.Sqrt2));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>float precision floating point number</returns>
        public float Entropy
        {
            get
            {
                return 0.5f + 0.5f * Maths.Log(2 * Maths.Pi * sigma * sigma) + mu;
            }
        }
        #endregion
    }
}
