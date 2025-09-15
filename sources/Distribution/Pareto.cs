using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Pareto distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Pareto_distribution
    /// </remarks>
    [Serializable]
    public class Pareto : IDistribution
    {
        #region Private data
        private float xm = 1;
        private float k = 1;
        #endregion

        #region Pareto components
        /// <summary>
        /// Initializes the Pareto distribution.
        /// </summary>
        public Pareto() { }
        /// <summary>
        /// Initializes the Pareto distribution.
        /// </summary>
        /// <param name="xm">Scale factor θ (0, +inf)</param>
        /// <param name="k">Parameter k (0, +inf)</param>
        public Pareto(float xm, float k)
        {
            Xm = xm; K = k;
        }
        /// <summary>
        /// Gets or sets the scale factor Xm (0, +inf).
        /// </summary>
        public float Xm
        {
            get
            {
                return this.xm;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.xm = value;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor k (0, +inf).
        /// </summary>
        public float K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value <= 0)
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
                return new RangeFloat(xm, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                if (k > 1)
                {
                    return (k * xm) / (k - 1);
                }
                return float.PositiveInfinity;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        /// <remarks>
        /// For k ≤ 2, the variance diverges to infinity.
        /// </remarks>
        public float Variance
        {
            get
            {
                if (k > 2)
                {
                    float kMinus1 = k - 1;
                    return (k * xm * xm) / (kMinus1 * kMinus1 * (k - 2));
                }
                return float.PositiveInfinity;
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                return new float[] { xm };
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return xm * Maths.Sqrt(2, k);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                if (k > 3)
                {
                    return 2 * (1 + k) / (k - 3) * Maths.Sqrt((k - 2) / k);
                }
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value. Requires k > 4 for a finite result.
        /// Formula: 6 * (k^3 + k^2 - 6k - 2) / (k * (k - 3) * (k - 4))
        /// </remarks>
        public float Excess
        {
            get
            {
                if (k <= 2)
                {
                    return float.NaN;
                }
                if (k <= 4)
                {
                    return float.PositiveInfinity;
                }

                float k2 = k * k;
                float k3 = k2 * k;
                return 6 * (k3 + k2 - 6 * k - 2) / (k * (k - 3) * (k - 4));
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < xm)
            {
                return 0;
            }
            return k * Maths.Pow(xm, k) / Maths.Pow(x, k + 1);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < xm)
            {
                return 0;
            }
            return 1.0f - Maths.Pow(xm / x, k);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return Maths.Log(xm / k) + 1f + 1f / k;
            }
        }
        #endregion
    }
}
