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
        private float mu = 0;
        private float beta = 1;
        #endregion

        #region Gumbel components
        /// <summary>
        /// Initializes the Gumbel distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ ∈ (-inf, +inf)</param>
        /// <param name="beta">Scale factor β ∈ (0, +inf)</param>
        public Gumbel(float mu, float beta)
        {
            Mu = mu; Beta = beta;
        }
        /// <summary>
        /// Gets or sets the shift factor μ ∈ (-inf, +inf).
        /// </summary>
        public float Mu
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
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return mu + beta * MathsF.Gamma;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return mu - beta * (float)Math.Log(Math.Log(2));
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return (float)((Math.PI * Math.PI) / 6.0) * beta * beta; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return mu;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return 1.14f;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                return 12.0f / 5;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get { return (float)Math.Log(beta) + MathsF.Gamma + 1; }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            float z = (x - mu) / beta;
            return (float)Math.Exp(-Math.Exp(-z));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            float z = (x - mu) / beta;
            return (1 / beta) * (float)Math.Exp(-(z + Math.Exp(-z)));
        }
        #endregion
    }
}
