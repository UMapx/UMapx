using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous complex Gaussian wavelet.
    /// </summary>
    [Serializable]
    public class ComplexGaussianWavelet : IComplexWavelet
    {
        #region Private data
        private int derivative;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous complex Gaussian wavelet.
        /// </summary>
        /// <param name="derivative">Derivative order [1, 8]</param>
        public ComplexGaussianWavelet(int derivative = 1)
        {
            Derivative = derivative;
        }
        /// <summary>
        /// Gets or sets the derivative order [1, 8].
        /// </summary>
        public int Derivative
        {
            get
            {
                return derivative;
            }
            set
            {
                if (value < 1 || value > 8)
                    throw new ArgumentException("Invalid argument value");

                derivative = value;
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
            float x2 = x * x;
            Complex32 f0 = Maths.Exp(-x2);
            Complex32 f1 = Maths.Exp(-Maths.I * x);
            Complex32 f2 = (f1 * f0) / Math.Pow(Math.Exp(-0.5) * Math.Pow(2, 0.5) * Math.Pow(Maths.Pi, 0.5), 0.5);
            Complex32 psi = 0;

            // Gaussian wavelet ('):
            switch (derivative)
            {
                case 1:
                    psi = f2 * (-Maths.I - 2 * x) * Math.Pow(2, 0.5);
                    break;

                case 2:
                    psi = 1.0 / 3 * f2 * (-3 + 4 * Maths.I * x + 4 * x2) * Math.Pow(6, 0.5);
                    break;

                case 3:
                    psi = 1.0 / 15 * f2 * (7 * Maths.I + 18 * x - 12 * Maths.I * x * x - 8 * Math.Pow(x, 3)) * Math.Pow(30, 0.5);
                    break;

                case 4:
                    psi = 1.0 / 105 * f2 * (25 - 56 * Maths.I * x - 72 * x * x + 32 * Maths.I * Math.Pow(x, 3) + 16 * Math.Pow(x, 4)) * Math.Pow(210, 0.5);
                    break;

                case 5:
                    psi = 1.0 / 315 * f2 * (-81 * Maths.I - 250 * x + 280 * Maths.I * x * x + 240 * Maths.Pow(x, 3) - 80 * Maths.I * Math.Pow(x, 4) - 32 * Math.Pow(x, 5)) * Math.Pow(210, 0.5);
                    break;

                case 6:
                    psi = 1.0 / 3465 * f2 * (-331 + 972 * Maths.I * x + 1500 * x * x - 1120 * Maths.I * Math.Pow(x, 3) - 720 * Math.Pow(x, 4) + 192 * Maths.I * Math.Pow(x, 5) + 64 * Math.Pow(x, 6)) * Math.Pow(2310, 0.5);
                    break;

                case 7:
                    psi = 1.0 / 45045 * f2 * (1303 * Maths.I + 4634 * x - 6804 * Maths.I * x * x - 7000 * Math.Pow(x, 3) + 3920 * Maths.I * Math.Pow(x, 4) + 2016 * Math.Pow(x, 5) - 448 * Maths.I * Math.Pow(x, 6) - 128 * Math.Pow(x, 7)) * Maths.Pow(30030, 0.5);
                    break;

                case 8:
                    psi = 1.0 / 45045 * f2 * (5937 - 20848 * Maths.I * x - 37072 * x * x + 36288 * Maths.I * Math.Pow(x, 3) + 28000 * Math.Pow(x, 4) - 12544 * Maths.I * Math.Pow(x, 5) - 5376 * Math.Pow(x, 6) + 1024 * Maths.I * Math.Pow(x, 7) + 256 * Math.Pow(x, 8)) * Math.Pow(20021, 0.5);
                    break;
            }

            return psi;
        }
        #endregion
    }
}
