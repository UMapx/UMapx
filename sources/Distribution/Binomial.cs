using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the binomial distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Binomial_distribution
    /// </remarks>
    [Serializable]
    public class Binomial : IDistribution
    {
        #region Private data
        private int n = 20;
        private float p = 0.5f;
        private float q = 0.5f;
        #endregion

        #region Binomial components
        /// <summary>
        /// Initializes the binomial distribution.
        /// </summary>
        public Binomial() { }
        /// <summary>
        /// Initializes the binomial distribution.
        /// </summary>
        /// <param name="n">Number of experiments (>0)</param>
        /// <param name="p">Probability of success [0, 1]</param>
        public Binomial(int n, float p)
        {
            N = n; P = p;
        }
        /// <summary>
        /// Gets or sets number of experiments.
        /// </summary>
        public int N
        {
            get
            {
                return n;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.n = value;
            }
        }
        /// <summary>
        /// Gets or sets probability of success [0, 1].
        /// </summary>
        public float P
        {
            get
            {
                return p;
            }
            set
            {
                if (value > 1 || value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.p = value;
                this.q = 1.0f - p;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(0, this.n);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return n * p;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return n * p * q;
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                float test = (n + 1) * p;

                if (test <= 0f)
                {
                    return new float[] { 0f };
                }

                if (test >= n + 1f)
                {
                    return new float[] { n };
                }

                float floor = Maths.Floor(test);
                float rounded = Maths.Round(test);

                if (Maths.Abs(test - rounded) < 1e-6f)
                {
                    if (rounded <= 0f)
                    {
                        return new float[] { 0f };
                    }

                    if (rounded >= n + 1f)
                    {
                        return new float[] { n };
                    }

                    return new float[] { rounded - 1f, rounded };
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
                // See e.g. https://en.wikipedia.org/wiki/Binomial_distribution#Median
                float median = Maths.Floor((n + 1) * p);
                return median > n ? n : median;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        /// <remarks>
        /// Skewness is undefined when all trials succeed or fail.
        /// </remarks>
        public float Skewness
        {
            get
            {
                if (p == 0f || p == 1f)
                    return float.NaN;
                return (q - p) / Maths.Sqrt(n * p * q);
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
                return (1 - 6 * p * q) / Variance;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < 0 || x > n)
            {
                return 0;
            }

            float a = Special.LogBinomial(n, x);
            float b = x == 0 ? 0 : x * Maths.Log(p);
            float c = n - x;
            float d = Maths.Log(1 - p);
            float log = a + b + c * d;

            return Maths.Exp(log);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < 0)
                return 0;
            if (x >= n)
                return 1;

            // Interpret x as the integer number of successes k.
            int k = (int)Maths.Floor(x);
            float a = n - k;
            float b = k + 1;
            return Special.BetaIncompleteRegularized(a, b, q);
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
