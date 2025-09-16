using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the negative binomial distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Negative_binomial_distribution
    /// </remarks>
    [Serializable]
    public class NegativeBinomial : IDistribution
    {
        #region Private data
        private int r = 1;
        private float p = 0.5f;
        private float q = 0.5f;
        #endregion

        #region Negative binomial components
        /// <summary>
        /// Initializes the negative binomial distribution with r = 1 and p = 0.5.
        /// </summary>
        public NegativeBinomial()
        {
        }
        /// <summary>
        /// Initializes the negative binomial distribution.
        /// </summary>
        /// <param name="r">Number of required successes (&gt; 0)</param>
        /// <param name="p">Success probability in each trial (0, 1]</param>
        public NegativeBinomial(int r, float p)
        {
            R = r;
            P = p;
        }
        /// <summary>
        /// Gets or sets the required number of successes (&gt; 0).
        /// </summary>
        public int R
        {
            get
            {
                return r;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                r = value;
            }
        }
        /// <summary>
        /// Gets or sets the success probability in each trial (0, 1].
        /// </summary>
        public float P
        {
            get
            {
                return p;
            }
            set
            {
                if (value <= 0f || value > 1f)
                    throw new ArgumentException("Invalid argument value");

                p = value;
                q = 1f - p;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(0f, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                if (q == 0f)
                    return 0f;

                return r * q / p;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (q == 0f)
                    return 0f;

                return r * q / (p * p);
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                if (q == 0f || r <= 1)
                {
                    return new float[] { 0f };
                }

                float value = (r - 1f) * q / p;
                float rounded = Maths.Round(value);
                if (Maths.Abs(value - rounded) < 1e-6f)
                {
                    float first = Maths.Max(0f, rounded - 1f);
                    float second = Maths.Max(0f, rounded);

                    if (Maths.Abs(first - second) < 1e-6f)
                    {
                        return new float[] { first };
                    }

                    return new float[] { first, second };
                }

                float floor = Maths.Floor(value);
                if (floor < 0f)
                {
                    floor = 0f;
                }

                return new float[] { floor };
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                if (q == 0f)
                    return float.NaN;

                return (2f - p) / Maths.Sqrt(r * q);
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        public float Excess
        {
            get
            {
                if (q == 0f)
                    return float.NaN;

                return (6f - 6f * p + p * p) / (r * q);
            }
        }
        /// <summary>
        /// Returns the value of the probability mass function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < 0f)
            {
                return 0f;
            }

            int k = (int)Maths.Floor(x);
            if (x != k)
            {
                return 0f;
            }

            if (q == 0f)
            {
                return k == 0 ? 1f : 0f;
            }

            float log = Special.LogBinomial(k + r - 1, r - 1) + k * Maths.Log(q) + r * Maths.Log(p);
            return Maths.Exp(log);
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < 0f)
            {
                return 0f;
            }

            if (q == 0f)
            {
                return 1f;
            }

            int k = (int)Maths.Floor(x);
            return Special.BetaIncompleteRegularized(r, k + 1, p);
        }
        /// <summary>
        /// Gets the entropy value.
        /// </summary>
        public float Entropy
        {
            get
            {
                return float.NaN;
            }
        }
        #endregion
    }
}
