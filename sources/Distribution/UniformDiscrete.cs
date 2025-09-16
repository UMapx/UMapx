using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the discrete uniform distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_uniform_distribution
    /// </remarks>
    [Serializable]
    public class UniformDiscrete : IDistribution
    {
        #region Private data
        private int a = 0;
        private int b = 1;
        private int n = 2;
        #endregion

        #region Uniform discrete components
        /// <summary>
        /// Initializes the discrete uniform distribution on the interval [0, 1].
        /// </summary>
        public UniformDiscrete()
        {
        }
        /// <summary>
        /// Initializes the discrete uniform distribution.
        /// </summary>
        /// <param name="a">Lower bound (inclusive).</param>
        /// <param name="b">Upper bound (inclusive).</param>
        public UniformDiscrete(int a, int b)
        {
            this.a = a;
            this.b = b;
            UpdateN();
        }
        /// <summary>
        /// Gets or sets the lower bound (inclusive).
        /// </summary>
        public int A
        {
            get
            {
                return a;
            }
            set
            {
                a = value;
                UpdateN();
            }
        }
        /// <summary>
        /// Gets or sets the upper bound (inclusive).
        /// </summary>
        public int B
        {
            get
            {
                return b;
            }
            set
            {
                b = value;
                UpdateN();
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(a, b);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return 0.5f * (a + b);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                float nn = n;
                return (nn * nn - 1f) / 12f;
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                if (n == 1)
                {
                    return new float[] { a };
                }

                return new float[] { float.NaN };
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return 0.5f * (a + b);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                if (n <= 1)
                    return float.NaN;

                return 0f;
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        public float Excess
        {
            get
            {
                if (n <= 1)
                    return float.NaN;

                float nn = n;
                return -6f * (nn * nn + 1f) / (5f * (nn * nn - 1f));
            }
        }
        /// <summary>
        /// Returns the value of the probability mass function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Function(float x)
        {
            int k = (int)Maths.Floor(x);
            if (x != k)
            {
                return 0f;
            }

            if (k < a || k > b)
            {
                return 0f;
            }

            return 1f / n;
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Distribution(float x)
        {
            if (x < a)
            {
                return 0f;
            }

            int k = (int)Maths.Floor(x);
            if (k >= b)
            {
                return 1f;
            }

            int count = k - a + 1;
            return count / (float)n;
        }
        /// <summary>
        /// Gets the entropy value.
        /// </summary>
        public float Entropy
        {
            get
            {
                return Maths.Log(n);
            }
        }
        /// <summary>
        /// Updates the number of points in the support.
        /// </summary>
        private void UpdateN()
        {
            if (b < a)
                throw new ArgumentException("Upper bound must not be less than lower bound");

            n = b - a + 1;
        }
        #endregion
    }
}
