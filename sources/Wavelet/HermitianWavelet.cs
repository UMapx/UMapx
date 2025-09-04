using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Hermitian wavelet.
    /// </summary>
    [Serializable]
    public class HermitianWavelet : IComplexWavelet
    {
        #region Private data
        private int derivative;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Hermitian wavelet.
        /// </summary>
        /// <param name="derivative">Derivative order [1, 3]</param>
        public HermitianWavelet(int derivative = 1)
        {
            Derivative = derivative;
        }
        /// <summary>
        /// Gets or sets the derivative order [1, 3].
        /// </summary>
        public int Derivative
        {
            get
            {
                return derivative;
            }
            set
            {
                if (value < 1 || value > 3)
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
            Complex32 psi = 0;
            Complex32 f0 = Math.Pow(Math.PI, -0.25) * Maths.Exp(-x2 / 2);

            switch (derivative)
            {
                case 1:
                    psi = Math.Sqrt(2) * x * f0;
                    break;

                case 2:
                    psi = 2.0 * Math.Sqrt(3.0) / 3.0 * f0 * (1 - x2);
                    break;

                case 3:
                    psi = 2.0 * Math.Sqrt(30.0) / 15.0 * f0 * (x2 * x - 3 * x);
                    break;
            }

            return psi;
        }
        #endregion
    }
}
