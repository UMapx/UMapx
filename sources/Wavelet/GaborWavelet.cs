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
        private double x0;
        private double k0;
        private double a;
        private double a2;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous complex Gabor wavelet.
        /// </summary>
        /// <param name="x0">Initial value</param>
        /// <param name="k0">Modulation factor</param>
        /// <param name="a">Factor</param>
        public GaborWavelet(double x0 = 0, double k0 = 1, double a = 2)
        {
            X0 = x0; K0 = k0; A = a;
        }
        /// <summary>
        /// Gets or sets the initial value.
        /// </summary>
        public double X0
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
        public double K0
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
        public double A
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
            double d = x - x0;
            return Math.Exp(-d * d / a2) * Maths.Exp(-Maths.I * k0 * d);
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public double WaveletReal(double x)
        {
            return Wavelet(x).Real;
        }
        #endregion
    }
}
