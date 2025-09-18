using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Burr distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Burr_distribution
    /// </remarks>
    [Serializable]
    public class Burr : IDistribution
    {
        #region Private data
        private float c;
        private float k;
        #endregion

        #region Burr distribution
        /// <summary>
        /// Initializes the Burr distribution.
        /// </summary>
        /// <param name="c">Form parameter c > 0</param>
        /// <param name="k">Scale parameter k > 0</param>
        public Burr(float c, float k)
        {
            C = c; K = k;
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter c > 0.
        /// </summary>
        public float C
        {
            get
            {
                return this.c;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.c = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter k > 0.
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
        /// Gets the mean value.
        /// Mean does not exist for k ≤ 1/c.
        /// </summary>
        public float Mean
        {
            get
            {
                if (k <= 1.0f / c)
                    return float.PositiveInfinity; // mean does not exist for k ≤ 1/c

                return k * Special.Beta(k - 1.0f / c, 1.0f + 1.0f / c);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        /// <remarks>
        /// The variance is defined only for <c>k &gt; 2/c</c>. For
        /// <c>k ≤ 2/c</c> the second moment diverges and the method
        /// returns <see cref="float.PositiveInfinity"/>.
        /// </remarks>
        public float Variance
        {
            get
            {
                if (k > 2.0f / c)
                {
                    float m1 = k * Special.Beta(k - 1.0f / c, 1.0f + 1.0f / c);
                    float m2 = k * Special.Beta(k - 2.0f / c, 1.0f + 2.0f / c);
                    return m2 - m1 * m1;
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
                if (c > 1)
                    return new float[] { Maths.Pow((c - 1) / (k * c + 1), 1.0f / c) };

                return new float[] { 0f };
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return Maths.Pow(Maths.Pow(2, 1.0f / k) - 1.0f, 1.0f / c);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        /// <remarks>
        /// Defined only when <c>k &gt; 3/c</c>; otherwise returns
        /// <see cref="float.NaN"/>.
        /// </remarks>
        public float Skewness
        {
            get
            {
                if (k > 3.0f / c)
                {
                    float m1 = k * Special.Beta(k - 1.0f / c, 1.0f + 1.0f / c);
                    float m2 = k * Special.Beta(k - 2.0f / c, 1.0f + 2.0f / c);
                    float m3 = k * Special.Beta(k - 3.0f / c, 1.0f + 3.0f / c);

                    float mu = m1;
                    float var = m2 - mu * mu;
                    float mu3 = m3 - 3.0f * mu * m2 + 2.0f * mu * mu * mu;

                    return mu3 / (var * Maths.Sqrt(var));
                }

                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value. The excess is defined only
        /// for <c>k &gt; 4/c</c>; for <c>2/c &lt; k ≤ 4/c</c> the fourth moment
        /// diverges and the method returns <see cref="float.PositiveInfinity"/>.
        /// For <c>k ≤ 2/c</c> the result is <see cref="float.NaN"/>.
        /// </remarks>
        public float Excess
        {
            get
            {
                if (k > 4.0f / c)
                {
                    float m1 = k * Special.Beta(k - 1.0f / c, 1.0f + 1.0f / c);
                    float m2 = k * Special.Beta(k - 2.0f / c, 1.0f + 2.0f / c);
                    float m3 = k * Special.Beta(k - 3.0f / c, 1.0f + 3.0f / c);
                    float m4 = k * Special.Beta(k - 4.0f / c, 1.0f + 4.0f / c);

                    float mu = m1;
                    float mu2 = m2 - mu * mu;
                    float mu4 = m4 - 4.0f * mu * m3 + 6.0f * mu * mu * m2 - 3.0f * mu * mu * mu * mu;

                    return mu4 / (mu2 * mu2) - 3.0f;
                }

                if (k > 2.0f / c)
                    return float.PositiveInfinity;

                return float.NaN;
            }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <remarks>
        /// Defined for all <c>c &gt; 0</c> and <c>k &gt; 0</c>.
        /// </remarks>
        public float Entropy
        {
            get
            {
                float dig1 = Special.DiGamma(1.0f);
                float digk = Special.DiGamma(k);

                return -Maths.Log(c * k) - (c - 1.0f) / c * (dig1 - digk) + (k + 1.0f) / k;
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
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= 0)
            {
                return 0;
            }

            return 1.0f - Maths.Pow(1.0f + Maths.Pow(x, c), -k);
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
                return 0;
            }

            float a = c * k;
            float b = Maths.Pow(x, c - 1);
            float d = 1 + Maths.Pow(x, c);
            return a * b / Maths.Pow(d, k + 1);
        }
        #endregion
    }
}
