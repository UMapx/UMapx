using System;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Haar wavelet.
    /// </summary>
    [Serializable]
    public class HaarWavelet : IDoubleWavelet
    {
        #region Haar wavelet components
        /// <summary>
        /// Initializes the continuous Haar wavelet.
        /// </summary>
        public HaarWavelet() { }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public double Scaling(double x)
        {
            if (0 <= x && x < 1)
            {
                return 1.0;
            }
            return 0.0;
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public double Wavelet(double x)
        {
            if (x >= 0)
            {
                return (x < 0.5) ? 1.0 : -1.0;
            }
            return 0.0;
        }
        #endregion
    }
}
