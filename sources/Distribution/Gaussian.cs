using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Gaussian distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Normal_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Gaussian : IDistribution
    {
        #region Private data
        private double sigma = 1;
        private double mu = 0;
        #endregion;

        #region Gaussian components
        /// <summary>
        /// Initializes the Gaussian distribution.
        /// </summary>
        public Gaussian() { }
        /// <summary>
        /// Initializes the Gaussian distribution.
        /// </summary>
        /// <param name="sigma">Standard deviation</param>
        /// <param name="mu">Mathematical expectation</param>
        public Gaussian(double sigma, double mu)
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
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Gets or sets the mathematical expectation.
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
                return this.mu;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return sigma * sigma;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return this.mu;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return this.mu;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0.0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 0.0;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return Maths.Exp(Maths.Pow((x - mu), 2) / (-2.0 * sigma * sigma)) / (Maths.Sqrt(2.0 * Maths.Pi) * sigma);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            return 0.5 + 0.5 * Special.Erf((x - mu) / Maths.Sqrt(2.0 * sigma * sigma));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(sigma * Maths.Sqrt(2 * Maths.Pi * Maths.E), Maths.E);
            }
        }
        #endregion
    }
}
