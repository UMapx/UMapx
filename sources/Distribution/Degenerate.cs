using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the degenerate distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Degenerate_distribution
    /// </remarks>
    [Serializable]
    public class Degenerate : IDistribution
    {
        #region Private data
        private int value;
        #endregion

        #region Degenerate components
        /// <summary>
        /// Initializes the degenerate distribution at zero.
        /// </summary>
        public Degenerate()
        {
        }
        /// <summary>
        /// Initializes the degenerate distribution.
        /// </summary>
        /// <param name="value">Single supported value.</param>
        public Degenerate(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Gets or sets the unique value whose probability is equal to one.
        /// </summary>
        public int Value
        {
            get
            {
                return this.value;
            }
            set
            {
                this.value = value;
            }
        }
        /// <summary>
        /// Gets the support interval consisting of a single point.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(value, value);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return value;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return 0f;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float[] Mode
        {
            get
            {
                return new float[] { value };
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return value;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        public float Excess
        {
            get
            {
                return float.NaN;
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
            return k == value ? 1f : 0f;
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Distribution(float x)
        {
            int k = (int)Maths.Floor(x);
            return k < value ? 0f : 1f;
        }
        /// <summary>
        /// Gets the entropy value.
        /// </summary>
        public float Entropy
        {
            get
            {
                return 0f;
            }
        }
        #endregion
    }
}
