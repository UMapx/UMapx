namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the discrete wavelet.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wavelet
    /// </remarks>
    public partial class WaveletPack
    {
        #region Haar wavelet
        /// <summary>
        /// Returns Haar wavelet.
        /// </summary>
        public static WaveletPack Haar
        {
            get
            {
                return Create(new float[] { 0.707106781186548f, 0.707106781186548f });
            }
        }
        #endregion
    }
}
