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
        #region Fbsp wavelets
        /// <summary>
        /// Returns B-spline wavelet 1-0-0.
        /// </summary>
        /// <remarks>
        /// Haar wavelet (delayed).
        /// </remarks>
        public static WaveletPack Fbsp100
        {
            get
            {
                return Create(new float[] { 0, 0, 0, 0, 0.707106781186548f, -0.707106781186548f, 0, 0, 0, 0 });
            }
        }
        /// <summary>
        /// Returns B-spline wavelet 1-0-3.
        /// </summary>
        public static WaveletPack Fbsp103
        {
            get
            {
                return Create(new float[] { -0.044194173824159f, 0.044194173824159f, 0.707106781186548f, 0.707106781186548f, 0.044194173824159f, -0.044194173824159f });
            }
        }
        /// <summary>
        /// Returns B-spline wavelet 1-0-5.
        /// </summary>
        public static WaveletPack Fbsp105
        {
            get
            {
                return Create(new float[] { 0.008286407592030f, -0.008286407592030f, -0.060766989008219f, 0.060766989008219f, 0.707106781186548f, 0.707106781186548f, 0.060766989008219f, -0.060766989008219f, -0.008286407592030f, 0.008286407592030f });
            }
        }
        #endregion
    }
}
