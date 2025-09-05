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
        #region Vaidyanathan wavelet
        /// <summary>
        /// Returns Vaidyanathan wavelet.
        /// </summary>
        public static WaveletPack Vaidyanathan
        {
            get
            {
                // Vaidyanathan wavelet:
                return WaveletPack.Create(new float[]
                {
                     -0.0000629061180000f, 0.0003436319050000f, -0.0004539566200000f,-0.0009448971360000f,0.0028438345470000f,0.0007081375040000f,-0.0088391034090000f,0.0031538470560000f,0.0196872150100000f,-0.0148534480050000f,-0.0354703986070000f,0.0387426192930000f,0.0558925236910000f,-0.0777097509020000f,-0.0839288843660000f,0.1319716614170000f,0.1350842271290000f,-0.1944504717660000f,-0.2634948024880000f,0.2016121617750000f,0.6356010598720000f,0.5727977932110000f,0.2501841295050000f,0.0457993341110000f
                });
            }
        }
        #endregion
    }
}
