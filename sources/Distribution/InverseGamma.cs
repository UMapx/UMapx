using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the inverse gamma distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Inverse-gamma_distribution
    /// </remarks>
    [Serializable]
    public class InverseGamma : IDistribution
    {
        #region Private data
        private float alpha = 1f;
        private float beta = 1f;
        private float constant = 1f;
        #endregion

        #region Inverse gamma components
        /// <summary>
        /// Initializes the inverse gamma distribution.
        /// </summary>
        public InverseGamma()
        {
            UpdateConstant();
        }
        /// <summary>
        /// Initializes the inverse gamma distribution.
        /// </summary>
        /// <param name="alpha">Shape parameter α (0, +inf).</param>
        /// <param name="beta">Scale parameter β (0, +inf).</param>
        public InverseGamma(float alpha, float beta)
        {
            Alpha = alpha;
            Beta = beta;
        }
        /// <summary>
        /// Gets or sets the shape parameter α (0, +inf).
        /// </summary>
        public float Alpha
        {
            get { return this.alpha; }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.alpha = value;
                UpdateConstant();
            }
        }
        /// <summary>
        /// Gets or sets the scale parameter β (0, +inf).
        /// </summary>
        public float Beta
        {
            get { return this.beta; }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.beta = value;
                UpdateConstant();
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(float.Epsilon, float.PositiveInfinity); }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                if (alpha <= 1f)
                    return float.PositiveInfinity;
                return beta / (alpha - 1f);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (alpha <= 2f)
                    return float.PositiveInfinity;
                float denom = (alpha - 1f) * (alpha - 1f) * (alpha - 2f);
                return (beta * beta) / denom;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { return float.NaN; }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get { return new float[] { beta / (alpha + 1f) }; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { return alpha > 3f ? (4f * Maths.Sqrt(alpha - 2f)) / (alpha - 3f) : float.NaN; }
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
                if (alpha <= 4f)
                    return float.NaN;

                float a = alpha;
                float num = 6f * (a * a * a + a * a - 6f * a - 2f);
                float den = a * (a - 3f) * (a - 4f);
                return num / den;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= 0f)
                return 0f;

            return constant * Maths.Pow(x, -alpha - 1f) * Maths.Exp(-beta / x);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= 0f)
                return 0f;

            return Special.GammaQ(alpha, beta / x);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get { return alpha + Maths.Log(beta) + Special.LogGamma(alpha) - (1f + alpha) * Special.DiGamma(alpha); }
        }
        #endregion

        #region Private methods
        private void UpdateConstant()
        {
            this.constant = Maths.Pow(beta, alpha) / Special.Gamma(alpha);
        }
        #endregion
    }
}
