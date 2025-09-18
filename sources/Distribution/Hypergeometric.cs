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

                return d * (k / n) * ((n - k) / n) * ((n - d) / (n - 1f));
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                float value = (k + 1f) * (d + 1f) / (n + 2f);
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
                double population = n;
                double successes = k;
                double draws = d;

                double numerator = population * population * (population - 1d) *
                                   (population * (population + 1d) - 6d * successes * (population - successes) -
                                    6d * draws * (population - draws)) +
                                   6d * draws * successes * (population - successes) * (population - draws) *
                                   (5d * population - 6d);
                double denominator = draws * successes * (population - successes) * (population - draws) *
                                     (population - 2d) * (population - 3d);

                if (denominator <= 0d)
                {
                    return numerator > 0d ? float.PositiveInfinity : float.NaN;
                }

                return (float)(numerator / denominator);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        /// <remarks>
        /// The function is defined only for integer values of <paramref name="x"/>.
        /// </remarks>
        public float Function(float x)
        {
            int k = (int)Math.Floor(x);
            if (x != k)
            {
                return 0f;
            }

            if (k < Math.Max(0, this.k + d - n) || k > Math.Min(d, this.k))
            {
                return 0f;
            }

            double a = Special.Binomial(d, k);
            double b = Special.Binomial(n - d, this.k - k);
            double c = Special.Binomial(n, this.k);
            return (float)(a * b / c);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            var support = Support;
            float lower = support.Min;
            float upper = support.Max;

            if (x < lower)
            {
                return 0f;
            }

            if (x >= upper)
            {
                return 1f;
            }

            float sum = 0f;
            int k = (int)Math.Floor(x);

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
