using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Meyer wavelet.
    /// </summary>
    [Serializable]
    public class MeyerWavelet : IFloatWavelet
    {
        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Meyer wavelet.
        /// </summary>
        public MeyerWavelet() { }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public float Scaling(float x)
        {
            // 2015, Victor Vermehren Valenzuela and H. M. de Oliveira gave 
            // the explicit expressions of Meyer wavelet and scale functions:
            if (x == 0)
            {
                return 2.0f / 3 + 4.0f / (3 * Maths.Pi);
            }
            float phiupper = Maths.Sin(2.0f * Maths.Pi / 3 * x) + 4.0f / 3 * x * Maths.Cos(4 * Maths.Pi / 3 * x);
            float phidown = Maths.Pi * x - 16 * Maths.Pi / 9 * Maths.Pow(x, 3);
            return phiupper / phidown;
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public float Wavelet(float x)
        {
            // 2015, Victor Vermehren Valenzuela and H. M. de Oliveira gave 
            // the explicit expressions of Meyer wavelet and scale functions:
            //
            // Kernel value:
            float t = x - 0.5f;
            // Finding ψ1(t):
            float psi1upper = 4.0f / (3 * Maths.Pi) * t * Maths.Cos(2 * Maths.Pi / 3 * t) - 1.0f / Maths.Pi * Maths.Sin(4 * Maths.Pi / 3 * t);
            float psi1down = t - 16.0f / 9 * Maths.Pow(t, 3);
            // Finding ψ2(t):
            float psi2upper = 8.0f / (3 * Maths.Pi) * t * Maths.Cos(8 * Maths.Pi / 3 * t) + 1.0f / Maths.Pi * Maths.Sin(4 * Maths.Pi / 3 * t);
            float psi2down = t - 64.0f / 9 * Maths.Pow(t, 3);
            // Finding ψ(t) = ψ1(t) + ψ2(t):
            return psi1upper / psi1down + psi2upper / psi2down;
        }
        #endregion
    }
}
