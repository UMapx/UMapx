using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Fisher distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/F-distribution
    /// </remarks>
    [Serializable]
    public class FisherSnedecor : IDistribution
    {
        #region Private data
        private int d1;
        private int d2;
        #endregion

        #region Fisher-Snedecor components
        /// <summary>
        /// Initializes the Fisher distribution.
        /// </summary>
        /// <param name="d1">First degree of freedom</param>
        /// <param name="d2">Second degree of freedom</param>
        public FisherSnedecor(int d1 = 1, int d2 = 1)
        {
            this.D1 = d1;
            this.D2 = d2;
        }
        /// <summary>
        /// Gets the value of the first degree of freedom.
        /// </summary>
        public int D1
        {
            get 
            {
                return d1; 
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentOutOfRangeException(nameof(d1), "The value must be greater than zero");

                d1 = value;
            }
        }
        /// <summary>
        /// Gets the value of the second degree of freedom.
        /// </summary>
        public int D2
        {
            get 
            { 
                return d2; 
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentOutOfRangeException(nameof(d2), "The value must be greater than zero");

                d2 = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        /// <remarks>
        /// The mean exists only for <c>d2 &gt; 2</c>; otherwise,
        /// <see cref="float.PositiveInfinity"/> is returned.
        /// </remarks>
        public float Mean
        {
            get
            {
                if (d2 <= 2)
                {
                    return float.PositiveInfinity;
                }

                return d2 / (d2 - 2.0f);
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
        /// <remarks>
        /// The variance is finite only when <c>d2 &gt; 4</c>.
        /// It becomes <see cref="float.PositiveInfinity"/> for
        /// <c>d2 &lt;= 4</c>; in particular, for <c>d2 &lt;= 2</c> the mean is
        /// already infinite.
        /// </remarks>
        public float Variance
        {
            get
            {
                if (d2 <= 2)
                {
                    // Mean is already infinite
                    return float.PositiveInfinity;
                }

                if (d2 <= 4)
                {
                    return float.PositiveInfinity;
                }

                return 2.0f * d2 * d2 * (d1 + d2 - 2) / (d1 * (d2 - 2) * (d2 - 2) * (d2 - 4));
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        /// <remarks>
        /// When <c>d1 ≤ 2</c>, the mode occurs at zero.
        /// </remarks>
        public float[] Mode
        {
            get
            {
                if (d1 > 2)
                {
                    float a = (d1 - 2.0f) / d1;
                    float b = d2 / (d2 + 2.0f);
                    return new float[] { a * b };
                }

                return new float[] { 0f };
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                if (d2 > 6)
                {
                    float v1 = 2 * d1 + d2 - 2;
                    float v2 = Maths.Sqrt(8 * (d2 - 4));
                    float v3 = Maths.Sqrt(d1 * (d1 + d2 - 2));
                    float v4 = d2 - 6;
                    return v1 * v2 / (v3 * v4);
                }
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value. The excess is defined
        /// only for <c>d2 &gt; 8</c>; otherwise, the fourth central moment
        /// diverges and <see cref="float.PositiveInfinity"/> is returned.
        /// </remarks>
        public float Excess
        {
            get
            {
                if (d2 <= 8)
                {
                    return float.PositiveInfinity;
                }

                float num = 12f * (d1 * (5f * d2 - 22f) * (d1 + d2 - 2f) + (d2 - 4f) * (d2 - 6f));
                float den = d1 * (d2 - 2f) * (d2 - 4f) * (d2 - 6f) * (d2 - 8f);
                return num / den;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(0, float.PositiveInfinity); }
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
            if (x <= 0)
                return 0;

            float u = (d1 * x) / (d1 * x + d2);
            return Special.BetaIncompleteRegularized(d1 * 0.5f, d2 * 0.5f, u);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= 0)
                return 0;
            
            float b = Special.Beta(d1 * 0.5f, d2 * 0.5f);
            float u = Maths.Pow(d1 * x, d1) * Maths.Pow(d2, d2) /
                Maths.Pow(d1 * x + d2, d1 + d2);
            return Maths.Sqrt(u) / (x * b);
        }
        #endregion
    }
}
