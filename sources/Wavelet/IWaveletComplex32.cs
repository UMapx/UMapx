using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the interface for continuous complex wavelets.
    /// </summary>
    public interface IWaveletComplex32
    {
        #region Interface
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        Complex32 Scaling(float x);
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        Complex32 Wavelet(float x);
        #endregion
    }
}
