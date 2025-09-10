using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Hermitian hat wavelet.
    /// </summary>
    [Serializable]
    public class HermitianHatWavelet : IWaveletComplex32
    {
        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Hermitian Hat wavelet.
        /// </summary>
        public HermitianHatWavelet() { }
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
            float cf = 2.0f / Maths.Sqrt(5) * Maths.Pow(Maths.Pi, -0.25f);
            return cf * (1 - x2 + Maths.I * x) * Maths.Exp(-0.5f * x2);
        }
        #endregion
    }
}
