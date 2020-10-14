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
        private double sigma = 1;
        private double mu = 0;
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
        public LogGaussian(double sigma, double mu)
        {
            Sigma = sigma;
            Mu = mu;
        }
        /// <summary>
        /// Gets or sets the standard deviation.
        /// </summary>
        public double Sigma
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
            get
            {
                return Maths.Exp(mu + sigma * sigma / 2);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return (Maths.Exp(sigma * sigma) - 1) * Maths.Exp(2 * mu + sigma * sigma);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Maths.Exp(mu);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return Maths.Exp(mu - sigma * sigma);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (Maths.Exp(sigma * sigma) + 2.0) * Maths.Sqrt(Maths.Exp(sigma * sigma) - 1.0);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return Maths.Exp(4 * sigma * sigma) + 2.0 * Maths.Exp(3 * sigma * sigma) + 3.0 * Maths.Exp(3 * sigma * sigma) - 6.0;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Maths.Exp(Maths.Pow((Maths.Log(x) - mu), 2) / (-2.0 * sigma * sigma)) / (Maths.Sqrt(2.0 * Maths.Pi) * sigma * x);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 0.5 + 0.5 * Special.Erf((Maths.Log(x) - mu) / Maths.Sqrt(sigma * 1.414));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return 0.5 + 0.5 * Maths.Log(2 * Maths.Pi * sigma * sigma) + mu;
            }
        }
        #endregion
    }
}
