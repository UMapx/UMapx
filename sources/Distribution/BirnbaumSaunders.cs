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

        #region Birnbaum-Saunders components
        /// <summary>
        /// Initializes the Birnbaum-Saunders distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ ∈ (-inf, +inf)</param>
        /// <param name="beta">Scale factor β ∈ (0, +inf)</param>
        /// <param name="gamma">Shape factor γ ∈ (0, +inf)</param>
        public BirnbaumSaunders(float mu, float beta, float gamma)
        {
            Mu = mu; Beta = beta; Gamma = gamma;
        }
        /// <summary>
        /// Gets or sets the shift factor μ ∈ (-inf, +inf).
        /// </summary>
        /// <example>
        /// <code>
        /// var distribution = new BirnbaumSaunders(-2f, 1f, 0.5f);
        /// float probability = distribution.Distribution(-1f); // returns 0.5
        /// </code>
        /// </example>
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
                return (float)((double)mu + beta);
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
        /// Gets the mode values.
        /// </summary>
        /// <remarks>
        /// Derived from solving a cubic equation; valid for γ &gt; 0.
        /// </remarks>
        public float[] Mode
        {
            get
            {
                double g2 = (double)gamma * gamma;
                double low = 0, high = g2 > 1 ? 1 / (3 * g2 - 1) : 1;
                for (int i = 0; i < 80; i++)
                {
                    double t = (low + high) / 2;
                    double value = ((t + 1 + g2) * t + 3 * g2 - 1) * t - 1;
                    if (value > 0) high = t; else low = t;
                }
                return new[] { (float)(mu + beta * ((low + high) / 2)) };
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        /// <remarks>
        /// See the standard moment expressions for the Birnbaum–Saunders distribution.
        /// </remarks>
        public float Skewness
        {
            get
            {
                float g2 = gamma * gamma;
                float numerator = 4f * gamma * (11f * g2 + 6f);
                float denominator = Maths.Pow(5f * g2 + 4f, 1.5f);
                return numerator / denominator;
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// See the standard moment expressions for the Birnbaum–Saunders distribution.
        /// </remarks>
        public float Excess
        {
            get
            {
                float g2 = gamma * gamma;
                float numerator = 6f * g2 * (93f * g2 + 40f);
                float denominator = 5f * g2 + 4f;
                denominator *= denominator;
                return numerator / denominator;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        /// <remarks>
        /// Computed from the normal-variable transformation for gamma &gt; 0.
        /// </remarks>
        public float Entropy
        {
            get
            {
                double g = gamma;
                double correction = DistributionNumerics.Integrate(z =>
                    Math.Exp(-0.5 * z * z) / Math.Sqrt(2 * Math.PI) *
                    DistributionNumerics.Log1p(g * g * z * z / 4), 0, 12);
                return (float)(0.5 + DistributionNumerics.LogSqrtTwoPi + Math.Log(beta) + Math.Log(g) - correction);
            }
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
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
