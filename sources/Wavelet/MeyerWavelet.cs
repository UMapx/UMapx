using System;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Meyer wavelet.
    /// </summary>
    [Serializable]
    public class MeyerWavelet : IDoubleWavelet
    {
        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Meyer wavelet.
        /// </summary>
        public MeyerWavelet() { }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public double Scaling(double x)
        {
            // 2015, Victor Vermehren Valenzuela and H. M. de Oliveira gave 
            // the explicit expressions of Meyer wavelet and scale functions:
            if (x == 0)
            {
                return 2.0 / 3 + 4.0 / (3 * Math.PI);
            }
            double phiupper = Math.Sin(2.0 * Math.PI / 3 * x) + 4.0 / 3 * x * Math.Cos(4 * Math.PI / 3 * x);
            double phidown = Math.PI * x - 16 * Math.PI / 9 * Math.Pow(x, 3);
            return phiupper / phidown;
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public double Wavelet(double x)
        {
            // 2015, Victor Vermehren Valenzuela and H. M. de Oliveira gave 
            // the explicit expressions of Meyer wavelet and scale functions:
            //
            // Kernel value:
            double t = x - 0.5;
            // Finding ψ1(t):
            double psi1upper = 4.0 / (3 * Math.PI) * t * Math.Cos(2 * Math.PI / 3 * t) - 1.0 / Math.PI * Math.Sin(4 * Math.PI / 3 * t);
            double psi1down = t - 16.0 / 9 * Math.Pow(t, 3);
            // Finding ψ2(t):
            double psi2upper = 8.0 / (3 * Math.PI) * t * Math.Cos(8 * Math.PI / 3 * t) + 1.0 / Math.PI * Math.Sin(4 * Math.PI / 3 * t);
            double psi2down = t - 64.0 / 9 * Math.Pow(t, 3);
            // Finding ψ(t) = ψ1(t) + ψ2(t):
            return psi1upper / psi1down + psi2upper / psi2down;
        }
        #endregion
    }
}
