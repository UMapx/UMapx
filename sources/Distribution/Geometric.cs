using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the geometric distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Geometric_distribution
    /// </remarks>
    [Serializable]
    public class Geometric : IDistribution
    {
        #region Private data
        private float p = 0.2f;
        private float q = 0.8f;
        #endregion

        #region Geometric components
        /// <summary>
        /// Initializes the geometric distribution.
        /// </summary>
        public Geometric() { }
        /// <summary>
        /// Initializes the geometric distribution.
        /// </summary>
        /// <param name="p">Probability of "success" (0, 1]</param>
        public Geometric(float p)
        {
            P = p;
        }
        /// <summary>
        /// Gets or sets the probability value of "success" (0, 1].
        /// </summary>
        public float P
        {
            get
            {
                return this.p;
            }
            set
            {
                if (value <= 0 || value > 1)
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
                return new RangeFloat(0, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return q / p;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return q / p / p;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return 0;
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
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return (2 - p) / Maths.Sqrt(1 - p);
            }
        }
        /// <summary>
        /// Returns the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get
            {
                return p * p / q + 3f;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < 0)
            {
                return 0f;
            }
            int k = (int)Maths.Floor(x);
            if (x != k)
            {
                return 0f;
            }
            return Maths.Pow(q, k) * p;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Number of failures before the first success</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < 0)
            {
                return 0;
            }
            x = Maths.Floor(x); // x is interpreted as the number of failures before first success
            return 1 - Maths.Pow(q, x + 1);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// When <c>p = 1</c>, the distribution is degenerate and the entropy equals zero.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                if (p == 1f)
                {
                    return 0f;
                }
                return -Maths.Log2(p) - q / p * Maths.Log2(q);
            }
        }
        #endregion
    }
}
