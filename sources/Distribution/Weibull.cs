using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Weibull distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Weibull_distribution
    /// </remarks>
    [Serializable]
    public class Weibull : IDistribution
    {
        #region Private data
        private float l = 1;
        private float k = 1;
        #endregion

        #region Weibull components
        /// <summary>
        /// Initializes the Weibull distribution.
        /// </summary>
        public Weibull() { }
        /// <summary>
        /// Initializes the Weibull distribution.
        /// </summary>
        /// <param name="lambda">Scale factor (0, + inf)</param>
        /// <param name="k">Shape factor (0, + inf)</param>
        public Weibull(float lambda, float k)
        {
            Lambda = lambda;
            K = k;
        }
        /// <summary>
        /// Gets or sets the value of the scale factor (0, + inf).
        /// </summary>
        public float Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.l = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the form factor (0, + inf).
        /// </summary>
        public float K
        {
            get
            {
                return this.k;
            }
            set
            {
                this.k = Maths.Max(0.0000001f, value);
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
                return l * Special.Gamma(1.0f + 1.0f / k);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return l * l * (Special.Gamma(1f + 2f / k) - Maths.Pow(Special.Gamma(1f + 1f / k), 2));
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return k > 1 ? l * Maths.Pow((k - 1f) / k, 1f / k) : 0f;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return l * Maths.Pow(Maths.Log(2), 1.0f / k);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return (Special.Gamma(1.0f + 3.0f / k) * Maths.Pow(l, 3) - 3 * Mean * Special.Gamma(1 + 2.0f / k) * Maths.Pow(l, 2) + 2 * Maths.Pow(Mean, 3)) / (Variance * Maths.Sqrt(Variance));
            }
        }
        /// <summary>
        /// Gets the excess coefficient (kurtosis minus 3).
        /// </summary>
        public float Excess
        {
            get
            {
                return (Special.Gamma(1.0f + 4.0f / k) * Maths.Pow(l, 4) - 4 * Mean * Special.Gamma(1 + 3.0f / k) * Maths.Pow(l, 3) + 6 * Maths.Pow(Mean, 2) * Maths.Pow(l, 2) * Special.Gamma(1.0f + 2.0f / k) - 3 * Maths.Pow(Mean, 4)) / (Variance * Variance) - 3f;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < 0)
            {
                return 0;
            }
            return (k / l) * Maths.Pow(x / l, k - 1) * Maths.Exp(-Maths.Pow(x / l, k));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 1 - Maths.Exp(-Maths.Pow(x / l, k));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return Maths.Gamma * (1f - 1f / k) + Maths.Log(l / k) + 1f;
            }
        }
        #endregion
    }
}
