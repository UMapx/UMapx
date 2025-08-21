namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the interface for continuous wavelets.
    /// </summary>
    public interface IFloatWavelet
    {
        #region Interface
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        float Scaling(float x);
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        float Wavelet(float x);
        #endregion
    }
}
