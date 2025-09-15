using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the uniform distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
    /// </remarks>
    [Serializable]
    public class Uniform : IDistribution
    {
        #region Private data
        private float a = 1;
        private float b = 1;
        #endregion

        #region Uniform components
        /// <summary>
        /// Initializes the uniform distribution.
        /// </summary>
        public Uniform() { }
        /// <summary>
        /// Initializes the uniform distribution.
        /// </summary>
        /// <param name="a">Shift parameter a</param>
        /// <param name="b">Shift parameter b</param>
        public Uniform(float a, float b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the shift parameter a.
        /// </summary>
        public float A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the shift parameter b.
        /// </summary>
        public float B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument [a, b].
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
                return (a + b) / 2.0f;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return Maths.Pow(b - a, 2) / 12.0f;
            }
        }
        /// <summary>
        /// Gets the mode values. Since the uniform distribution has no unique mode,
        /// the value is undefined and represented by <see cref="float.NaN"/>.
        /// </summary>
        public float[] Mode
        {
            get
            {
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
                return Mean;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return 0;
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
                return -1.2f;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < a)
            {
                return 0;
            }
            else if (x > b)
            {
                return 0;
            }
            return 1.0f / (b - a);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < a)
            {
                return 0;
            }
            else if (x > b)
            {
                return 1;
            }
            return (x - a) / (b - a);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return Maths.Log(b - a);
            }
        }
        #endregion
    }
}
