using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Burr distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Burr_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Burr : IDistribution
    {
        #region Private data
        private float c;
        private float k;
        #endregion

        #region Burr distribution
        /// <summary>
        /// Initializes the Burr distribution.
        /// </summary>
        /// <param name="c">Form parameter c > 0</param>
        /// <param name="k">Scale parameter k > 0</param>
        public Burr(float c, float k)
        {
            C = c; K = k;
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter c > 0.
        /// </summary>
        public float C
        {
            get
            {
                return this.c;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.c = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter k > 0.
        /// </summary>
        public float K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return k * Special.Beta(k - 1.0f / c, 1.0f + 1.0f / c); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return (float)Math.Pow((c - 1) / (k * c + 1), 1.0 / c);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return (float)Math.Pow(Math.Pow(2, 1.0 / k) - 1.0, 1.0 / c);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(float.Epsilon, float.PositiveInfinity); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= 0)
            {
                return float.NaN;
            }

            return 1.0f - (float)Math.Pow(1.0 + Math.Pow(x, c), -k);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= 0)
            {
                return float.NaN;
            }

            float a = c * k;
            float b = (float)Math.Pow(x, c - 1);
            float d = 1 + (float)Math.Pow(x, c);
            return a * b / (float)Math.Pow(d, k + 1);
        }
        #endregion
    }
}
