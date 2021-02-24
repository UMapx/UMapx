using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous complex Morlet wavelet.
    /// </summary>
    [Serializable]
    public class ComplexMorletWavelet : IComplexWavelet
    {
        #region Private data
        private float fb;
        private float fc;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous complex Morlet wavelet.
        /// </summary>
        /// <param name="fb">Bandwidth</param>
        /// <param name="fc">Center frequency</param>
        public ComplexMorletWavelet(float fb = 0.5f, float fc = 1)
        {
            Fb = fb; Fc = fc;
        }
        /// <summary>
        /// Gets or sets the bandwidth.
        /// </summary>
        public float Fb
        {
            get
            {
                return this.fb;
            }
            set
            {
                this.fb = value;
            }
        }
        /// <summary>
        /// Gets or sets the center frequency.
        /// </summary>
        public float Fc
        {
            get
            {
                return this.fc;
            }
            set
            {
                this.fc = value;
            }
        }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public Complex32 Scaling(float x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public Complex32 Wavelet(float x)
        {
            return Math.Pow(Maths.Pi * fb, -0.5) * Maths.Exp(2 * Maths.Pi * Maths.I * fc * x) * Math.Exp(-(x * x) / fb);
        }
        #endregion
    }
}
