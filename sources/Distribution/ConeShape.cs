using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the distribution of the conical shape.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Cone-shape_distribution_function
    /// </remarks>
    [Serializable]
    public class ConeShape : IDistribution
    {
        #region Private data
        private float a;
        #endregion

        #region Cone-Shape components
        /// <summary>
        /// Initializes the distribution of the conical shape.
        /// </summary>
        /// <param name="a">Coefficient</param>
        public ConeShape(float a = 0.001f)
        {
            A = a;
        }
        /// <summary>
        /// Gets or sets the coefficient value.
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
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get { throw new NotSupportedException(); }
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
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the kernel density function.
        /// </summary>
        /// <param name="eta">Value</param>
        /// <param name="tau">Value</param>
        /// <returns>Value</returns>
        public float Function(float eta, float tau)
        {
            float ksi = Maths.Pi * eta * tau;
            float psi = Maths.Exp(-2 * Maths.Pi * a * tau * tau);
            return Special.Sinc(ksi, 1) * psi;
        }
        /// <summary>
        /// Returns the value of the kernel distribution function.
        /// </summary>
        /// <param name="t">Value</param>
        /// <param name="tau">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float t, float tau)
        {
            if (tau == 0)
            {
                return t == 0 ? float.PositiveInfinity : 0f;
            }
            if (Math.Abs(tau) >= 2 * Math.Abs(t))
            {
                return 1.0f / tau * Maths.Exp(-2 * Maths.Pi * a * tau * tau);
            }
            return 0;
        }
        #endregion
    }
}
