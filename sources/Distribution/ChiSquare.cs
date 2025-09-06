using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the chi-square distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Chi-squared_distribution
    /// </remarks>
    [Serializable]
    public class ChiSquare : IDistribution
    {
        #region Private data
        private int k = 1;
        #endregion

        #region Chi-square components
        /// <summary>
        /// Initializes the chi-square distribution.
        /// </summary>
        /// <param name="k">Degrees of freedom (0, +inf)</param>
        public ChiSquare(int k)
        {
            K = k;
        }
        /// <summary>
        /// Gets or sets the degrees of freedom (0, +inf).
        /// </summary>
        public int K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(0, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return k;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return k - 2.0f / 3;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                if (k >= 2)
                {
                    return k - 2;
                }
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return 2 * k;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return (float)Math.Sqrt(8.0 / k);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                return 12.0f / k;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get
            {
                float k2 = k / 2.0f;
                float s1 = (float)Math.Log(2.0 * Special.Gamma(k2));
                float s2 = (1.0f - k2) * Special.DiGamma(k2);
                return k2 + s1 + s2;
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= 0)
            {
                return 0;
            }
            return Special.GammaP(k / 2.0f, x / 2.0f);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= 0)
            {
                return 1;
            }
            return Special.GammaQ(k / 2.0f, x / 2.0f);
        }
        #endregion
    }
}
