namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the wavelet transform interface.
    /// </summary>
    public interface IWaveletTransform
    {
        #region Interface
        /// <summary>
        /// Gets or sets the discrete wavelet decomposition.
        /// </summary>
        WaveletDecomposition WaveletDecomposition { get; set; }
        #endregion
    }
}
