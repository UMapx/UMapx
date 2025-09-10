using System;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Haar wavelet on [0,1).
    /// </summary>
    [Serializable]
    public class HaarWavelet : IWaveletFloat
    {
        #region Haar wavelet components
        /// <summary>
        /// Initializes the continuous Haar wavelet.
        /// </summary>
        public HaarWavelet() { }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Value</param>
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
        /// Returns the value of the wavelet function defined on [0,1).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public float Wavelet(float x)
        {
            if (0 <= x && x < 0.5f) return 1f;
            if (0.5f <= x && x < 1f) return -1f;
            return 0f;
        }
        #endregion
    }
}
