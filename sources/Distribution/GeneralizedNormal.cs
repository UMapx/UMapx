using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the generalized normal distribution (version 1 / exponential power distribution).
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Generalized_normal_distribution
    /// </remarks>
    [Serializable]
    public class GeneralizedNormal : IDistribution
    {
        #region Private data
        private float mu = 0f;
        private float alpha = 1f;
        private float beta = 2f;
        #endregion

        #region Generalized normal components
        /// <summary>
        /// Initializes the generalized normal distribution.
        /// </summary>
        public GeneralizedNormal() { }
        /// <summary>
        /// Initializes the generalized normal distribution.
        /// </summary>
        /// <param name="mu">Location parameter.</param>
        /// <param name="alpha">Scale parameter (0, +inf).</param>
        /// <param name="beta">Shape parameter (0, +inf).</param>
        public GeneralizedNormal(float mu, float alpha, float beta)
        {
            Mu = mu;
            Alpha = alpha;
            Beta = beta;
        }
        /// <summary>
        /// Gets or sets the location parameter.
        /// </summary>
        public float Mu
        {
            get { return this.mu; }
            set { this.mu = value; }
        }
        /// <summary>
        /// Gets or sets the scale parameter (0, +inf).
        /// </summary>
        public float Alpha
        {
            get { return this.alpha; }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.alpha = value;
            }
        }
        /// <summary>
        /// Gets or sets the shape parameter (0, +inf).
        /// </summary>
        public float Beta
        {
            get { return this.beta; }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.beta = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity); }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return this.mu; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                float ratio = Special.Gamma(3f / beta) / Special.Gamma(1f / beta);
                return alpha * alpha * ratio;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { return this.mu; }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get { return new float[] { this.mu }; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { return 0f; }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get
            {
                float g1 = Special.Gamma(5f / beta);
                float g2 = Special.Gamma(1f / beta);
                float g3 = Special.Gamma(3f / beta);
                float kurtosis = (g1 * g2) / (g3 * g3);
                return kurtosis - 3f;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            float z = Maths.Abs((x - mu) / alpha);
            float norm = beta / (2f * alpha * Special.Gamma(1f / beta));
            return norm * Maths.Exp(-Maths.Pow(z, beta));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            float z = x - mu;
            if (z == 0f)
                return 0.5f;

            float sign = z > 0 ? 1f : -1f;
            float u = Maths.Pow(Maths.Abs(z) / alpha, beta);
            float p = Special.GammaP(1f / beta, u);
            return 0.5f + 0.5f * sign * p;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                float a = 1f / beta;
                float norm = beta / (2f * alpha * Special.Gamma(a));
                return a - Maths.Log(norm);
            }
        }
        #endregion
    }
}
