using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the binomial distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// <see href="https://en.wikipedia.org/wiki/Binomial_distribution"/>.
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
        /// <param name="n">Number of experiments (nonnegative).</param>
        /// <param name="p">Probability of success [0, 1].</param>
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
        /// <remarks>Returns the lower median, rounded to binary32.</remarks>
        public float Median
        {
            get
            {
                if (p == 0 || n == 0) return 0;
                if (p == 1) return n;
                // Symmetry gives an exact lower median, including odd trial counts.
                if (p == 0.5f) return n / 2;
                int low = 0, high = n;
                while (low < high)
                {
                    int mid = low + (high - low) / 2;
                    double cumulative = Special.DistributionBeta(n - (double)mid, mid + 1.0, 1.0 - p);
                    if (cumulative >= 0.5) high = mid; else low = mid + 1;
                }
                return low;
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
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Function(float x)
        {
            if (float.IsNaN(x)) return float.NaN;
            if (x < 0 || (double)x > n || x != Math.Floor(x)) return 0;
            if (n == 0 || p == 0) return x == 0 ? 1 : 0;
            if (p == 1) return (double)x == n ? 1 : 0;
            double k = x;
            double log = Special.DistributionLogGamma(n + 1.0) - Special.DistributionLogGamma(k + 1)
                - Special.DistributionLogGamma(n - k + 1) + k * Math.Log(p)
                + (n - k) * DistributionNumerics.Log1p(-(double)p);
            return (float)Math.Exp(log);
        }
        /// <summary>
        /// Returns the value of the probability mass cumulative function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Distribution(float x)
        {
            if (float.IsNaN(x)) return float.NaN;
            if (x < 0) return 0;
            if ((double)x >= n) return 1;
            if (p == 0) return 1;
            if (p == 1) return 0;
            double k = Math.Floor(x);
            return (float)Special.DistributionBeta(n - k, k + 1, 1.0 - p);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value.</returns>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
}
