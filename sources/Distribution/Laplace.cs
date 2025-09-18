using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Laplace distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Laplace_distribution
    /// </remarks>
    [Serializable]
    public class Laplace : IDistribution
    {
        #region Private data
        private float a = 1;
        private float b = 0;
        #endregion

        #region Laplace components
        /// <summary>
        /// Initializes the Laplace distribution.
        /// </summary>
        public Laplace() { }
        /// <summary>
        /// Initializes the Laplace distribution.
        /// </summary>
        /// <param name="alfa">Scale factor (0, + inf)</param>
        /// <param name="beta">Shift coefficient</param>
        public Laplace(float alfa, float beta)
        {
            Alfa = alfa;
            Beta = beta;
        }
        /// <summary>
        /// Gets or sets the value of the scale factor (0, + inf).
        /// </summary>
        public float Alfa
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the shift coefficient.
        /// </summary>
        public float Beta
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
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return b;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return 2.0f / a / a;
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                return new float[] { b };
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return b;
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
                return 3;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            return a / 2.0f * Maths.Exp(-a * Maths.Abs(x - b));
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= b)
            {
                return 0.5f * Maths.Exp(a * (x - b));
            }
            return 1 - 0.5f * Maths.Exp(-a * (x - b));
        }
        /// <summary>
        /// Returns the value of differential entropy H = 1 - ln(a / 2).
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return 1 - Maths.Log(a / 2.0f);
            }
        }
        #endregion
    }
}
