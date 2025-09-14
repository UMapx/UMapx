using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Burr distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Burr_distribution
    /// </remarks>
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
        /// Mean does not exist for k ≤ 1/c.
        /// </summary>
        public float Mean
        {
            get
            {
                if (k <= 1.0f / c)
                    return float.PositiveInfinity; // mean does not exist for k ≤ 1/c

                return k * Special.Beta(k - 1.0f / c, 1.0f + 1.0f / c);
            }
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
                if (c > 1)
                    return Maths.Pow((c - 1) / (k * c + 1), 1.0f / c);

                return 0;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return Maths.Pow(Maths.Pow(2, 1.0f / k) - 1.0f, 1.0f / c);
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
            get { return new RangeFloat(0, float.PositiveInfinity); }
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
                return 0;
            }

            return 1.0f - Maths.Pow(1.0f + Maths.Pow(x, c), -k);
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
                return 0;
            }

            float a = c * k;
            float b = Maths.Pow(x, c - 1);
            float d = 1 + Maths.Pow(x, c);
            return a * b / Maths.Pow(d, k + 1);
        }
        #endregion
    }
}
