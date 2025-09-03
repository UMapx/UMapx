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
        private float mu = 0;
        private float beta = 1;
        private float gamma = 1;
        #endregion

        #region Birnbaum-Saunders components
        /// <summary>
        /// Initializes the Birnbaum-Saunders distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ ∈ (0, +inf)</param>
        /// <param name="beta">Scale factor β ∈ (0, +inf)</param>
        /// <param name="gamma">Shape factor γ ∈ (0, +inf)</param>
        public BirnbaumSaunders(float mu, float beta, float gamma)
        {
            Mu = mu; Beta = beta; Gamma = gamma;
        }
        /// <summary>
        /// Gets or sets the shift factor μ ∈ (0, +inf).
        /// </summary>
        public float Mu
        {
            get
            {
                return mu;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor β ∈ (0, +inf).
        /// </summary>
        public float Beta
        {
            get
            {
                return beta;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.beta = value;
            }
        }
        /// <summary>
        /// Gets or sets the form factor γ ∈ (0, +inf).
        /// </summary>
        public float Gamma
        {
            get
            {
                return gamma;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.gamma = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(this.mu, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return 1 + 0.5f * gamma * gamma;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return gamma * gamma * (1 + (5 * gamma * gamma) / 4); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get { throw new NotSupportedException(); }
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
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            float a = (float)Math.Sqrt(x);
            float b = (float)Math.Sqrt(1.0 / x);
            float z = (a - b) / gamma;

            // Normal cumulative distribution function
            return Special.Erfc(-z / 1.4142135623731f) * 0.5f;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            float c = x - mu;

            float a = (float)Math.Sqrt(c / beta);
            float b = (float)Math.Sqrt(beta / c);

            float alpha = (a + b) / (2 * gamma * c);
            float z = (a - b) / gamma;

            // Normal cumulative distribution function
            return alpha * Special.Erfc(-z / 1.4142135623731f) * 0.5f;
        }
        #endregion
    }
}
