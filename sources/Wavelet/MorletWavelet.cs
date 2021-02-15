using System;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Morlet wavelet.
    /// </summary>
    [Serializable]
    public class MorletWavelet : IFloatWavelet
    {
        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Morlet wavelet.
        /// </summary>
        public MorletWavelet()
        {
        }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public float Scaling(float x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public float Wavelet(float x)
        {
            float x2 = x * x;
            return (float)Math.Exp(-x2 / 2) * (float)Math.Cos(5 * x);
        }
        #endregion
    }
}
