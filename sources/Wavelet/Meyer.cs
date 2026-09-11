using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Meyer wavelet.
    /// </summary>
    [Serializable]
    public class Meyer : IWaveletFloat
    {
        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Meyer wavelet.
        /// </summary>
        public Meyer() { }
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
            if (Math.Abs(x) == 0.75f) return (float)(2 / (3 * Math.PI));
            // Double intermediates avoid cancellation at float samples adjacent
            // to the removable singularities; exact limits handle the poles.
            double t = x;
            double phiupper = Math.Sin(2 * Math.PI / 3 * t) + 4.0 / 3 * t * Math.Cos(4 * Math.PI / 3 * t);
            double phidown = Math.PI * t * (1 - 16.0 / 9 * t * t);
            return (float)(phiupper / phidown);
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
            double t = (double)x - 0.5;
            if (t == 0) return (float)(4 / Math.PI);
            // Finding ψ1(t):
            double psi1upper = 4.0 / (3 * Math.PI) * t * Math.Cos(2 * Math.PI / 3 * t) - Math.Sin(4 * Math.PI / 3 * t) / Math.PI;
            double psi1down = t * (1 - 16.0 / 9 * t * t);
            // Finding ψ2(t):
            double psi2upper = 8.0 / (3 * Math.PI) * t * Math.Cos(8 * Math.PI / 3 * t) + Math.Sin(4 * Math.PI / 3 * t) / Math.PI;
            double psi2down = t * (1 - 64.0 / 9 * t * t);
            // Finding ψ(t) = ψ1(t) + ψ2(t):
            double psi1 = Math.Abs(t) == 0.75 ? -1.0 / 3 : psi1upper / psi1down;
            double psi2 = Math.Abs(t) == 0.375 ? 4 / (3 * Math.PI) : psi2upper / psi2down;
            return (float)(psi1 + psi2);
        }
        #endregion
    }
}
