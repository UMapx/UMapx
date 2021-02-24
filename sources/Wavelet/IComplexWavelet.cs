using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the interface for continuous complex wavelets.
    /// </summary>
    public interface IComplexWavelet
    {
        #region Interface
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        Complex32 Scaling(float x);
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        Complex32 Wavelet(float x);
        #endregion
    }
}
