using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Gaussian wavelet.
    /// </summary>
    [Serializable]
    public class Gaussian : IWaveletFloat
    {
        #region Private data
        private int derivative;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Gaussian wavelet.
        /// </summary>
        /// <param name="derivative">Derivative order [1, 8]</param>
        public Gaussian(int derivative = 1)
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
        public float Scaling(float x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public float Wavelet(float x)
        {
            float x2 = x * x;
            float f0 = Maths.Pow(2.0f / Maths.Pi, 0.25f) * Maths.Exp(-x2);
            float psi = 0;

            // Gaussian wavelet ('):
            switch (derivative)
            {
                case 1:
                    psi = -2.0f * x * f0;
                    break;

                case 2:
                    psi = 2.0f / Maths.Pow(3, 0.5f) * (-1.0f + 2 * x2) * f0;
                    break;

                case 3:
                    psi = 4.0f / Maths.Pow(15, 0.5f) * x * (3 - 2 * x2) * f0;
                    break;

                case 4:
                    psi = 4.0f / Maths.Pow(105, 0.5f) * (3 - 12 * x2 + 4 * x2 * x2) * f0;
                    break;

                case 5:
                    psi = 8.0f / (3 * Maths.Pow(105, 0.5f)) * x * (-15 + 20 * x2 - 4 * x2 * x2) * f0;
                    break;

                case 6:
                    psi = 8.0f / (3 * Maths.Pow(1155, 0.5f)) * (-15 + 90 * x2 - 60 * x2 * x2 + 8 * Maths.Pow(x2, 3)) * f0;
                    break;

                case 7:
                    psi = 16.0f / (3 * Maths.Pow(15015, 0.5f)) * x * (105 - 210 * x2 + 84 * x2 * x2 - 8 * Maths.Pow(x2, 3)) * f0;
                    break;

                case 8:
                    psi = 16.0f / (45 * Maths.Pow(1001, 0.5f)) * (105 - 840 * x2 + 840 * x2 * x2 - 224 * Maths.Pow(x2, 3) + 16 * Maths.Pow(x2, 4)) * f0;
                    break;
            }

            return psi;
        }
        #endregion
    }
}
