using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the distribution of Bernoulli.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Bernoulli_distribution
    /// </remarks>
    [Serializable]
    public class Bernoulli : IDistribution
    {
        #region Private data
        private float p = 0.2f;
        private float q = 0.8f;
        #endregion

        #region Bernoully components
        /// <summary>
        /// Initializes a Bernoulli distribution.
        /// </summary>
        public Bernoulli() { }
        /// <summary>
        /// Initializes a Bernoulli distribution.
        /// </summary>
        /// <param name="p">Probability of success [0, 1]</param>
        public Bernoulli(float p)
        {
            P = p;
        }
        /// <summary>
        /// Gets or sets the probability of success [0, 1].
        /// </summary>
        public float P
        {
            get
            {
                return this.p;
            }
            set
            {
                if (value < 0 || value > 1)
                    throw new ArgumentException("Invalid argument value");

                this.p = value;
                this.q = 1 - this.p;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(0, 1);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return p;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return p * q;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                if (q > p)
                {
                    return 0;
                }
                else if (q < p)
                {
                    return 1;
                }
                return 0.5f;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                if (q > p)
                {
                    return 0;
                }
                else if (q < p)
                {
                    return 1;
                }
                return 0.5f;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return (p - q) / Maths.Sqrt(p * q);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                return (6 * p * p - 6 * p + 1) / (p * (1 - p));
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x == 0)
            {
                return q;
            }
            return p;
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
            else if (x >= 1)
            {
                return 1;
            }
            return q;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return -q * Maths.Log(q) - p * Maths.Log(p);
            }
        }
        #endregion
    }
}
