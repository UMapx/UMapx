using System;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Haar wavelet.
    /// </summary>
    [Serializable]
    public class HaarWavelet : IFloatWavelet
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
        public float Scaling(float x)
        {
            if (0 <= x && x < 1)
            {
                return 1.0f;
            }
            return 0.0f;
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public float Wavelet(float x)
        {
            if (x >= 0)
            {
                return (x < 0.5) ? 1.0f : -1.0f;
            }
            return 0.0f;
        }
        #endregion
    }
}
