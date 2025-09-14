using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Birnbaum-Saunders distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Birnbaum–Saunders_distribution
    /// </remarks>
    [Serializable]
    public class BirnbaumSaunders : IDistribution
    {
        #region Private data
        private float mu = 0;
        private float beta = 1;
        private float gamma = 1;
        #endregion

        #region Utility
        /// <summary>
        /// Returns the real cube root of the given value.
        /// </summary>
        private static float Cbrt(float x)
        {
            return (x >= 0f) ? Maths.Pow(x, 1f / 3f) : -Maths.Pow(-x, 1f / 3f);
        }
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
                if (value <= 0)
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
                if (value <= 0)
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
                if (value <= 0)
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
                return mu + beta * (1 + 0.5f * gamma * gamma);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        /// <remarks>
        /// The expression is valid for γ &gt; 0.
        /// </remarks>
        public float Median
        {
            get
            {
                float g2 = gamma * gamma;
                float t = 0.5f * gamma + Maths.Sqrt(0.25f * g2 + 1f);
                return mu + beta * t * t;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return beta * beta * gamma * gamma * (1 + 5f * gamma * gamma / 4f); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        /// <remarks>
        /// Derived from solving a cubic equation; valid for γ &gt; 0.
        /// </remarks>
        public float Mode
        {
            get
            {
                float g2 = gamma * gamma;
                float a = 1f + g2;
                float b = 3f * g2 - 1f;
                float p = b - a * a / 3f;
                float q = 2f * a * a * a / 27f - a * b / 3f + 1f;
                float d = q * q / 4f + p * p * p / 27f;
                if (d < 0) return float.NaN;
                float sqrt = Maths.Sqrt(d);
                float u = Cbrt(-q / 2f + sqrt);
                float v = Cbrt(-q / 2f - sqrt);
                float t = u + v - a / 3f;
                return mu + beta * t;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        /// <remarks>
        /// The skewness exists for all γ &gt; 0.
        /// </remarks>
        public float Skewness
        {
            get
            {
                float g2 = gamma * gamma;
                float den = Maths.Pow(5f * g2 + 4f, 1.5f);
                return 4f * gamma * (11f * g2 + 6f) / den;
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        /// <remarks>
        /// The expression is valid for γ &gt; 0.
        /// </remarks>
        public float Excess
        {
            get
            {
                float g2 = gamma * gamma;
                float den = (5f * g2 + 4f);
                den *= den;
                return 6f * g2 * (93f * g2 + 40f) / den;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        /// <remarks>
        /// Closed-form expression for γ &gt; 0.
        /// </remarks>
        public float Entropy
        {
            get
            {
                float g2 = gamma * gamma;
                return 0.5f * (1f + Maths.Log(2f * Maths.Pi)) +
                       Maths.Log(beta * gamma / 2f) + 0.5f * g2;
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= mu)
                return 0;

            float z = (Maths.Sqrt((x - mu) / beta) - Maths.Sqrt(beta / (x - mu))) / gamma;

            return 0.5f * Special.Erfc(-z / Maths.Sqrt2);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= mu)
                return 0;

            float a = Maths.Sqrt((x - mu) / beta);
            float b = Maths.Sqrt(beta / (x - mu));
            float z = (a - b) / gamma;
            float phi = Maths.Exp(-0.5f * z * z) / Maths.Sqrt(2 * Maths.Pi);

            return (a + b) / (2 * gamma * (x - mu)) * phi;
        }
        #endregion
    }
}
