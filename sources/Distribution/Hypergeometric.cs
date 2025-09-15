using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the hypergeometric distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hypergeometric_distribution
    /// </remarks>
    [Serializable]
    public class Hypergeometric : IDistribution
    {
        #region Private data
        private float n = 30;
        private float k = 20;
        private float d = 20;
        #endregion

        #region Hypergeometric components
        /// <summary>
        /// Initializes the hypergeometric distribution.
        /// </summary>
        public Hypergeometric() { }
        /// <summary>
        /// Initializes the hypergeometric distribution.
        /// </summary>
        /// <param name="n">Parameter N (0, +inf]</param>
        /// <param name="k">Parameter K [0, N]</param>
        /// <param name="d">Parameter D [0, N]</param>
        public Hypergeometric(float n, float k, float d)
        {
            N = n; K = k; D = d;
        }
        /// <summary>
        /// Gets or sets the value of the parameter N (0, +inf].
        /// </summary>
        public float N
        {
            get
            {
                return n;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.n = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter D [0, N].
        /// </summary>
        public float D
        {
            get
            {
                return d;
            }
            set
            {
                if (value < 0 || value > N)
                    throw new ArgumentException("Invalid argument value");

                this.d = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter K [0, N].
        /// </summary>
        public float K
        {
            get
            {
                return k;
            }
            set
            {
                if (value < 0 || value > N)
                    throw new ArgumentException("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument [max(0, K + D - N), min(D, K)].
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(Math.Max(0, k + d - n), Math.Min(d, k));
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                if (n <= 0f)
                    return 0f;

                return k * d / n;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (n <= 1f)
                    return 0f;

                return k * (d / n) * ((n - d) / n) * ((n - k) / (n - 1.0f));
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                float value = ((k + 1f) * (d + 1f)) / (n + 2f);
                float lower = Maths.Max(0f, k + d - n);
                float upper = Maths.Min(d, k);

                float floor = Maths.Range(Maths.Floor(value), lower, upper);
                float rounded = Maths.Round(value);

                if (Maths.Abs(value - rounded) < 1e-6f)
                {
                    float candidate1 = rounded - 1f;
                    float candidate2 = rounded;
                    bool firstValid = candidate1 >= lower - 1e-6f && candidate1 <= upper + 1e-6f;
                    bool secondValid = candidate2 >= lower - 1e-6f && candidate2 <= upper + 1e-6f;

                    if (firstValid && secondValid)
                    {
                        float mode1 = Maths.Range(candidate1, lower, upper);
                        float mode2 = Maths.Range(candidate2, lower, upper);

                        if (Maths.Abs(mode1 - mode2) < 1e-6f)
                        {
                            return new float[] { mode1 };
                        }

                        return new float[] { mode1, mode2 };
                    }

                    if (firstValid)
                    {
                        return new float[] { Maths.Range(candidate1, lower, upper) };
                    }

                    if (secondValid)
                    {
                        return new float[] { Maths.Range(candidate2, lower, upper) };
                    }
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
                float k1 = (n - 2 * d) * Maths.Pow(n - 1, 0.5f) * (n - 2 * k);
                float k2 = Maths.Sqrt(k * d * (n - d) * (n - k)) * (n - 2);
                return k1 / k2;
            }
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
                float n2 = n * n;
                float k1 = (n2 * (n - 1)) / (k * (n - 2) * (n - 3) * (n - k));
                float k2 = (n * (n + 1) - 6 * n * (n - k)) / (d * (n - d));
                float k3 = (3 * n * (n - k) * (n + 6)) / n2 - 6;
                return k1 * (k2 + k3);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < Math.Max(0, k + d - n) || x > Math.Min(d, k))
            {
                return 0;
            }

            double a = Special.Binomial(d, x);
            double b = Special.Binomial(n - d, k - x);
            double c = Special.Binomial(n, k);
            return (float)(a * b / c);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            float sum = 0;
            int k = (int)x;
            for (int i = 0; i <= k; i++)
            {
                sum += Function(i);
            }

            return sum;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
}
