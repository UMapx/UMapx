using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the beta distribution of the second kind.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Beta_prime_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class BetaPrime : IDistribution
    {
        #region Private data
        private double alpha = 1; // shape (α)
        private double beta = 1;  // shape (β)
        #endregion

        #region Beta-prime components
        /// <summary>
        /// Initializes beta distribution of the second kind.
        /// </summary>
        /// <param name="alpha">Parameter α (0, +inf)</param>
        /// <param name="beta">Parameter β (0, +inf)</param>
        public BetaPrime(double alpha, double beta)
        {
            Alpha = alpha; Beta = beta;
        }
        /// <summary>
        /// Gets or sets the value of the parameter α ∈ (0, +inf).
        /// </summary>
        public double Alpha
        {
            get
            {
                return alpha;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.alpha = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter β ∈ (0, +inf).
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
                if (beta > 1)
                {
                    return alpha / (beta - 1);
                }

                return Double.PositiveInfinity;
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
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (alpha >= 1)
                {
                    return (alpha - 1) / (beta + 1);
                }

                return 0.0;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                if (beta > 2.0)
                {
                    double num = alpha * (alpha + beta - 1);
                    double den = (beta - 2) * Math.Pow(beta - 1, 2);
                    return num / den;
                }
                else if (beta > 1.0)
                {
                    return Double.PositiveInfinity;
                }

                return double.NaN;
            }
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
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x <= 0)
            {
                return 0;
            }
            return Special.BetaIncomplete(alpha, beta, x / (1 + x));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x <= 0)
            {
                return 0;
            }

            double num = Math.Pow(x, alpha - 1) * Math.Pow(1 + x, -alpha - beta);
            double den = Special.Beta(alpha, beta);
            return num / den;
        }
        #endregion
    }
}
