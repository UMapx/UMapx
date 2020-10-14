using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Birnbaum-Saunders distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Birnbaum–Saunders_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class BirnbaumSaunders : IDistribution
    {
        #region Private data
        private double mu = 0;
        private double beta = 1;
        private double gamma = 1;
        #endregion

        #region Birnbaum-Saunders components
        /// <summary>
        /// Initializes the Birnbaum-Saunders distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ ∈ (0, +inf)</param>
        /// <param name="beta">Scale factor β ∈ (0, +inf).</param>
        /// <param name="gamma">Shape factor γ ∈ (0, +inf)</param>
        public BirnbaumSaunders(double mu, double beta, double gamma)
        {
            Mu = mu; Beta = beta; Gamma = gamma;
        }
        /// <summary>
        /// Gets or sets the shift factor μ ∈ (0, +inf).
        /// </summary>
        public double Mu
        {
            get
            {
                return mu;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

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
        /// Gets or sets the form factor γ ∈ (0, +inf).
        /// </summary>
        public double Gamma
        {
            get
            {
                return gamma;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.gamma = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(this.mu, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return 1 + 0.5 * gamma * gamma;
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
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return gamma * gamma * (1 + (5 * gamma * gamma) / 4); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
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
            double a = Math.Sqrt(x);
            double b = Math.Sqrt(1.0 / x);
            double z = (a - b) / gamma;

            // Normal cumulative distribution function
            return Special.Erfc(-z / 1.4142135623731) * 0.5;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double c = x - mu;

            double a = Math.Sqrt(c / beta);
            double b = Math.Sqrt(beta / c);

            double alpha = (a + b) / (2 * gamma * c);
            double z = (a - b) / gamma;

            // Normal cumulative distribution function
            return alpha * Special.Erfc(-z / 1.4142135623731) * 0.5;
        }
        #endregion
    }
}
