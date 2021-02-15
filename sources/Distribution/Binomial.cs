using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the binomial distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Binomial_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Binomial : IDistribution
    {
        #region Private data
        private float n = 20;
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
        public Binomial(float n, float p)
        {
            N = n; P = p;
        }
        /// <summary>
        /// Gets or sets number of experiments.
        /// </summary>
        public float N
        {
            get
            {
                return n;
            }
            set
            {
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
                    throw new Exception("Invalid argument value");

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
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                float test = (n + 1) * p;

                if (test <= 0 || (int)test != test)
                    return (float)Math.Floor(test);

                if (test <= n)
                    return test;

                if (test == n + 1)
                    return n;

                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return (float)Math.Floor(n * p);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return (q - p) / (float)Math.Sqrt(n * p * q);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
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
        /// <returns>float precision floating point number</returns>
        public float Function(float x)
        {
            if (x < 0 || x > n)
            {
                return 0;
            }

            float a = Special.LogBinomial(n, x);
            float b = x == 0 ? 0 : x * (float)Math.Log(p);
            float c = (n - x);
            float d = (float)Math.Log(1 - p);
            float log = a + b + c * d;

            return (float)Math.Exp(log);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>float precision floating point number</returns>
        public float Distribution(float x)
        {
            if (x < 0)
                return 0;
            if (x >= n)
                return 1;

            float a = n - x;
            float b = x + 1;
            return Special.Beta(a, b, q);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>float precision floating point number</returns>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
}
