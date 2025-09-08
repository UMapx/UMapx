using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the normalized continuous Morlet wavelet with zero mean.
    /// </summary>
    [Serializable]
    public class MorletWavelet : IFloatWavelet
    {
        #region Private data
        private float omega0 = 5;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Morlet wavelet.
        /// </summary>
        public MorletWavelet()
        {
        }
        /// <summary>
        /// Gets or sets the central frequency.
        /// </summary>
        public float Omega0
        {
            get
            {
                return this.omega0;
            }
            set
            {
                this.omega0 = value;
            }
        }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public float Scaling(float x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the normalized wavelet function with zero mean.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public float Wavelet(float x)
        {
            float x2 = x * x;
            return Maths.Exp(-x2 / 2) * (Maths.Cos(omega0 * x) - Maths.Exp(-0.5f * omega0 * omega0)) / Maths.Pow(Maths.Pi, 0.25f);
        }
        #endregion
    }
}
