using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the discrete wavelet.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wavelet
    /// </remarks>
    public partial class WaveletPacket
    {
        #region Kravchenko wavelet
        /// <summary>
        /// Returns Kravchenko wavelet.
        /// </summary>
        public static WaveletPacket Kravchenko
        {
            get
            {
                // Left values of scale function:
                float[] left = new float[]
                {
                    0.438708321041f,
                   -0.047099287129f,
                   -0.118027008279f,
                    0.037706980974f,
                    0.043603935723f,
                   -0.025214528289f,
                   -0.011459893503f,
                    0.013002207742f,
                   -0.001878954975f,
                   -0.003758906625f,
                    0.005085949920f,
                   -0.001349824585f,
                   -0.003639380570f,
                    0.002763059895f,
                    0.001188712844f,
                   -0.001940226446f,
                    0.000384982816f,
                    0.000499860951f,
                   -0.000700388155f,
                    0.000468702885f,
                    0.000255769244f,
                   -0.000649033581f,
                    0.000266223602f,
                    0.000307507863f,
                   -0.000463771747f,
                    0.000104807634f,
                    0.000324973138f,
                   -0.000288500372f,
                   -0.000066833177f,
                    0.000021430184f,
                   -0.000018524173f,
                   -0.000032851429f,
                   -0.000000000000f
                };

                // Right values of scale function:
                float[] right = new float[]
                {
                    0.757698251288f,
                    0.438708321041f,
                   -0.047099287129f,
                   -0.118027008279f,
                    0.037706980974f,
                    0.043603935723f,
                   -0.025214528289f,
                   -0.011459893503f,
                    0.013002207742f,
                   -0.001878954975f,
                   -0.003758906625f,
                    0.005085949920f,
                   -0.001349824585f,
                   -0.003639380570f,
                    0.002763059895f,
                    0.001188712844f,
                   -0.001940226446f,
                    0.000384982816f,
                    0.000499860951f,
                   -0.000700388155f,
                    0.000468702885f,
                    0.000255769244f,
                   -0.000649033581f,
                    0.000266223602f,
                    0.000307507863f,
                   -0.000463771747f,
                    0.000104807634f,
                    0.000324973138f,
                   -0.000288500372f,
                   -0.000066833177f,
                    0.000021430184f,
                   -0.000018524173f,
                   -0.000032851429f
                };

                // Kravchenko orthogonal wavelet:
                return WaveletPacket.Create(Matrice.Concat(left.Flip(), right));
            }
        }
        #endregion
    }
}
