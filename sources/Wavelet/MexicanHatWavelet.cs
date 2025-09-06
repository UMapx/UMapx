using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Mexican hat wavelet.
    /// </summary>
    [Serializable]
    public class MexicanHatWavelet : IFloatWavelet
    {
        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Mexican hat wavelet.
        /// </summary>
        public MexicanHatWavelet() { }
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
            return 2.0f / (float)(Math.Sqrt(3) * Maths.Pow(Maths.Pi, 0.25f)) * (1 - x2) * Maths.Exp(-x2 / 2);
        }
        #endregion
    }
}
