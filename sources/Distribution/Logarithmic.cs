using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the logarithmic distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Logarithmic_distribution
    /// </remarks>
    [Serializable]
    public class Logarithmic : IDistribution
    {
        #region Private data
        private float p = 0.66f;
        #endregion

        #region Logarithmic components
        /// <summary>
        /// Initializes the logarithmic distribution.
        /// </summary>
        /// <param name="p">Parameter (0 &lt; p &lt; 1)</param>
        public Logarithmic(float p)
        {
            P = p;
        }
        /// <summary>
        /// Gets or sets the value of the parameter p ∈ (0, 1).
        /// </summary>
        public float P
        {
            get
            {
                return p;
            }
            set
            {
                if (value <= 0 || value >= 1)
                    throw new ArgumentException("Invalid argument value");

                this.p = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(1, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return (float)(-1 / Math.Log(1 - p)) * p / (1 - p);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                float k1 = p + Maths.Log(1 - p);
                float k2 = Maths.Pow(1 - p, 2) * Maths.Pow(Maths.Log(1 - p), 2);
                return -p * k1 / k2;
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
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get { return new float[] { 1f }; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value (x ≥ 1)</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < 1)
            {
                return 0;
            }

            var n = Maths.Floor(x);
            float c = -1 / Maths.Log(1 - p);
            float sum = 0;

            for (int k = 1; k <= n; k++)
            {
                sum += c * Maths.Pow(p, k) / k;
            }

            return sum;
        }
        /// <summary>
        /// Returns the value of the probability mass function.
        /// </summary>
        /// <param name="x">Value (integer x ≥ 1)</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < 1)
            {
                return 0;
            }

            var k = Maths.Floor(x);

            if (x != k)
            {
                return 0;
            }

            return -1 / Maths.Log(1 - p) * Maths.Pow(p, k) / k;
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
}
