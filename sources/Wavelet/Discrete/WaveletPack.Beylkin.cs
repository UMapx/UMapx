namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the discrete wavelet.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wavelet
    /// </remarks>
    /// </summary>
    public partial class WaveletPack
    {
        #region Beylkin wavelet
        /// <summary>
        /// Returns Beylkin wavelet.
        /// </summary>
        public static WaveletPack Beylkin
        {
            get
            {
                // Beylkin wavelet:
                return WaveletPack.Create(new float[]
                {
                    0.0993057653740000f,0.4242153608130000f,0.6998252140570000f,0.4497182511490000f,-0.1109275983480000f,-0.2644972314460000f,0.0269003088040000f,0.1555387318770000f,-0.0175207462670000f,-0.0885436306230000f,0.0196798660440000f,0.0429163872740000f,-0.0174604086960000f,-0.0143658079690000f,0.0100404118450000f,0.0014842347820000f,-0.0027360316260000f,0.0006404853290000f
                });
            }
        }
        #endregion
    }
}
