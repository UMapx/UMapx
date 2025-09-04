using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Weibull distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Weibull_distribution
    /// </remarks>
    /// </summary>
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
                this.k = MathsF.Max(0.0000001f, value);
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
                return l * l * (Special.Gamma(1.0f + 2.0f / k) - Special.Gamma(1.0f + 1.0f / k));
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return l * MathsF.Pow(k - 1, 1.0f / k) / MathsF.Pow(k, 1.0f / k);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return l * MathsF.Pow(MathsF.Log(2), 1.0f / k);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return (Special.Gamma(1.0f + 3.0f / k) * MathsF.Pow(l, 3) - 3 * Mean * Special.Gamma(1 + 2.0f / k) * MathsF.Pow(l, 2) + 2 * MathsF.Pow(Mean, 3)) / (Variance * MathsF.Sqrt(Variance));
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                return (Special.Gamma(1.0f + 4.0f / k) * MathsF.Pow(l, 4) - 4 * Mean * Special.Gamma(1 + 3.0f / k) * MathsF.Pow(l, 3) + 6 * MathsF.Pow(Mean, 2) * MathsF.Pow(l, 2) * Special.Gamma(1.0f + 2.0f / k) - 3 * MathsF.Pow(Mean, 4)) / (Variance * Variance);
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
            return (k / l) * MathsF.Pow(x / l, k - 1) * MathsF.Exp(-MathsF.Pow(x / l, k));
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
            return 1 - MathsF.Exp(-MathsF.Pow(x / l, k));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return MathsF.Gamma * (1.0f - 1.0f / k) + MathsF.Pow(l / k, k) + MathsF.Log(l / k);
            }
        }
        #endregion
    }
}
