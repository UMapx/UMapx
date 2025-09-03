using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous complex Gabor wavelet.
    /// </summary>
    [Serializable]
    public class GaborWavelet : IComplexWavelet
    {
        #region Private data
        private float x0;
        private float k0;
        private float a;
        private float a2;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous complex Gabor wavelet.
        /// </summary>
        /// <param name="x0">Initial value</param>
        /// <param name="k0">Modulation factor</param>
        /// <param name="a">Factor</param>
        public GaborWavelet(float x0 = 0, float k0 = 1, float a = 2)
        {
            X0 = x0; K0 = k0; A = a;
        }
        /// <summary>
        /// Gets or sets the initial value.
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
        /// Gets or sets the modulation factor.
        /// </summary>
        public float K0
        {
            get
            {
                return this.k0;
            }
            set
            {
                this.k0 = value;
            }
        }
        /// <summary>
        /// Gets or sets the factor.
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
                this.a2 = a * a;
            }
        }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public Complex32 Scaling(float x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public Complex32 Wavelet(float x)
        {
            float d = x - x0;
            return Math.Exp(-d * d / a2) * MathF.Exp(-MathF.I * k0 * d);
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public float WaveletReal(float x)
        {
            return Wavelet(x).Real;
        }
        #endregion
    }
}
