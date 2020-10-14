namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the wavelet transform interface.
    /// </summary>
    public interface IWaveletTransform
    {
        #region Interface
        /// <summary>
        /// Gets or sets the discrete wavelet.
        /// </summary>
        WaveletPack Wavelet { get; set; }
        #endregion
    }
}
