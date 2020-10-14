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
        private double fb;
        private double fc;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous complex Morlet wavelet.
        /// </summary>
        /// <param name="fb">Bandwidth</param>
        /// <param name="fc">Center frequency</param>
        public ComplexMorletWavelet(double fb = 0.5, double fc = 1)
        {
            Fb = fb; Fc = fc;
        }
        /// <summary>
        /// Gets or sets the bandwidth.
        /// </summary>
        public double Fb
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
        public double Fc
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
        public Complex Scaling(double x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public Complex Wavelet(double x)
        {
            return Math.Pow(Maths.Pi * fb, -0.5) * Maths.Exp(2 * Maths.Pi * Maths.I * fc * x) * Math.Exp(-(x * x) / fb);
        }
        #endregion
    }
}
