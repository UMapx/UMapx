using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Wigner semicircular distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
    /// </remarks>
    [Serializable]
    public class Wigner : IDistribution
    {
        #region Private data
        private float r;
        #endregion

        #region Wigner components
        /// <summary>
        /// Initializes the Wigner semicircular distribution.
        /// </summary>
        /// <param name="r">Radius</param>
        public Wigner(float r)
        {
            R = r;
        }
        /// <summary>
        /// Gets or sets the radius value.
        /// </summary>
        public float R
        {
            get
            {
                return this.r;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.r = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(-r, r);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return r * r / 4.0f;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                return new float[] { 0f };
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
                return -1;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= -r || x >= r) return 0f;

            float r2 = r * r, x2 = x * x;
            float a = Maths.Sqrt(r2 - x2);
            float b = 2 / (float)(Math.PI * r2);
            return b * a;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= -r) return 0f;
            if (x >= r) return 1f;

            float r2 = r * r, x2 = x * x;
            float a = Maths.Sqrt(r2 - x2);
            float b = x / (Maths.Pi * r2);
            float c = Maths.Asin(x / r) / Maths.Pi;
            return 0.5f + b * a + c;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                // Standard differential entropy: H = 1/2 + ln(π r / 2)
                return Maths.Log(Maths.Pi * r / 2) + 0.5f;
            }
        }
        #endregion
    }
}
