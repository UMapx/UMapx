using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Cauchy distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Cauchy_distribution
    /// </remarks>
    [Serializable]
    public class Cauchy : IDistribution
    {
        #region Private data
        private float g = 0.5f;
        private float x0 = 0;
        #endregion

        #region Caushi components
        /// <summary>
        /// Initializes the Cauchy distribution.
        /// </summary>
        public Cauchy() { }
        /// <summary>
        /// Initializes the Cauchy distribution.
        /// </summary>
        /// <param name="gamma">Scale factor (0, + inf)</param>
        /// <param name="x0">Shift coefficient</param>
        public Cauchy(float gamma, float x0)
        {
            Gamma = gamma;
            X0 = x0;
        }
        /// <summary>
        /// Gets or sets the value of the scale factor (0, + inf).
        /// </summary>
        public float Gamma
        {
            get
            {
                return this.g;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.g = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the shift coefficient.
        /// </summary>
        public float X0
        {
            get
            {
                return this.x0;
            }
            set
            {
                this.x0 = value;
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
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return float.PositiveInfinity;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return x0;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return x0;
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
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            return 1.0f / (Maths.Pi * g * (1.0f + Maths.Pow((x - x0) / g)));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            return 1.0f / Maths.Pi * Maths.Atan((x - x0) / g) + 0.5f;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return Maths.Log(4 * Maths.Pi * g);
            }
        }
        #endregion
    }
}
