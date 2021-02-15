using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Wiener semicircular distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Wigner : IDistribution
    {
        #region Private data
        private float r;
        #endregion;

        #region Wigner components
        /// <summary>
        /// Initializes the Wiener semicircular distribution.
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
                    throw new Exception("Invalid argument value");

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
        /// Gets the kurtosis coefficient.
        /// </summary>
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
        /// <returns>float precision floating point number</returns>
        public float Function(float x)
        {
            if (Math.Abs(x) > r)
            {
                return float.NaN;
            }

            float r2 = r * r, x2 = x * x;
            float a = (float)Math.Sqrt(r2 - x2);
            float b = 2 / (float)(Math.PI * r2);
            return b * a;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>float precision floating point number</returns>
        public float Distribution(float x)
        {
            if (Math.Abs(x) > r)
            {
                return float.NaN;
            }

            float r2 = r * r, x2 = x * x;
            float a = (float)Math.Sqrt(r2 - x2);
            float b = x / (Maths.Pi * r2);
            float c = (float)Math.Asin(x / r) / Maths.Pi;
            return 0.5f + b * a + c;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>float precision floating point number</returns>
        public float Entropy
        {
            get
            {
                return Maths.Log(Maths.Pi * r) - 1.0f / 2;
            }
        }
        #endregion
    }
}
