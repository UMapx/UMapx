using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Gumbel distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gumbel_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Gumbel : IDistribution
    {
        #region Private data
        private double mu = 0;
        private double beta = 1;
        #endregion

        #region Gumbel components
        /// <summary>
        /// Initializes the Gumbel distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ ∈ (-inf, +inf)</param>
        /// <param name="beta">Scale factor β ∈ (0, +inf).</param>
        public Gumbel(double mu, double beta)
        {
            Mu = mu; Beta = beta;
        }
        /// <summary>
        /// Gets or sets the shift factor μ ∈ (-inf, +inf).
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
        /// Gets or sets the scale factor β ∈ (0, +inf).
        /// </summary>
        public double Beta
        {
            get
            {
                return beta;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.beta = value;
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
                return mu + beta * Maths.Gamma;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return mu - beta * Math.Log(Math.Log(2));
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return ((Math.PI * Math.PI) / 6.0) * beta * beta; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return mu;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 1.14;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 12.0 / 5;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { return Math.Log(beta) + Maths.Gamma + 1; }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            double z = (x - mu) / beta;
            return Math.Exp(-Math.Exp(-z));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double z = (x - mu) / beta;
            return (1 / beta) * Math.Exp(-(z + Math.Exp(-z)));
        }
        #endregion
    }
}
