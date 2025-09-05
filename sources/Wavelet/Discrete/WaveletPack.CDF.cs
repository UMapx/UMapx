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
        #region Cohen-Daubechies-Feaveau wavelets
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 1/1).
        /// </summary>
        public static WaveletPack CDF11
        {
            get
            {
                return WaveletPack.Bior11;
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 1/3).
        /// </summary>
        public static WaveletPack CDF13
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.707106781186548f,
                   0.707106781186548f }, new float[] {
                   0.088388347648318f,
                   0.088388347648318f,
                  -0.707106781186548f,
                   0.707106781186548f,
                  -0.088388347648318f,
                  -0.088388347648318f,
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 1/5).
        /// </summary>
        public static WaveletPack CDF15
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.707106781186548f,
                   0.707106781186548f }, new float[] {
                  -0.016572815184060f,
                  -0.016572815184060f,
                   0.121533978016438f,
                   0.121533978016438f,
                  -0.707106781186548f,
                   0.707106781186548f,
                  -0.121533978016438f,
                  -0.121533978016438f,
                   0.016572815184060f,
                   0.016572815184060f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 3/1).
        /// </summary>
        public static WaveletPack CDF31
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.176776695296637f,
                   0.530330085889911f,
                   0.530330085889911f,
                   0.176776695296637f }, new float[] {
                   0.353553390593274f,
                   1.060660171779821f,
                  -1.060660171779821f,
                  -0.353553390593274f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 5/1).
        /// </summary>
        public static WaveletPack CDF51
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.044194173824159f,
                   0.220970869120796f,
                   0.441941738241592f,
                   0.441941738241592f,
                   0.220970869120796f,
                   0.044194173824159f }, new float[] {
                  -0.265165042944955f,
                  -1.325825214724777f,
                  -1.767766952966369f,
                   1.767766952966369f,
                   1.325825214724777f,
                   0.265165042944955f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 5/3).
        /// </summary>
        public static WaveletPack CDF53
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.044194173824159f,
                   0.220970869120796f,
                   0.441941738241592f,
                   0.441941738241592f,
                   0.220970869120796f,
                   0.044194173824159f }, new float[] {
                  -0.055242717280199f,
                  -0.276213586400995f,
                  -0.817592215746946f,
                  -1.878252387526767f,
                  -2.099223256647564f,
                   1.436310649285175f,
                   0.773398041922786f,
                  -0.287262129857035f,
                  -0.276213586400995f,
                  -0.055242717280199f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 5/5).
        /// </summary>
        public static WaveletPack CDF55
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.044194173824159f,
                   0.220970869120796f,
                   0.441941738241592f,
                   0.441941738241592f,
                   0.220970869120796f,
                   0.044194173824159f }, new float[] {
                  -0.012084344405044f,
                  -0.060421722025218f,
                  -0.041432037960149f,
                   0.276213586400995f,
                   0.468527295932688f,
                  -0.543795498226959f,
                  -1.450121328605225f,
                   1.450121328605224f,
                   0.543795498226959f,
                  -0.468527295932688f,
                  -0.276213586400995f,
                   0.041432037960149f,
                   0.060421722025218f,
                   0.012084344405044f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 2/2).
        /// </summary>
        public static WaveletPack CDF22
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.353553390593274f,
                   0.707106781186548f,
                   0.353553390593274f,
                   }, new float[] {
                   0.176776695296637f,
                   0.353553390593274f,
                  -1.060660171779821f,
                   0.353553390593274f,
                   0.176776695296637f,
                   0.000000000000000f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 2/4).
        /// </summary>
        public static WaveletPack CDF24
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.353553390593274f,
                   0.707106781186548f,
                   0.353553390593274f,
                   }, new float[] {
                  -0.033145630368119f,
                  -0.066291260736239f,
                   0.176776695296637f,
                   0.419844651329513f,
                  -0.994368911043582f,
                   0.419844651329513f,
                   0.176776695296637f,
                  -0.066291260736239f,
                  -0.033145630368119f,
                   0.000000000000000f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 2/6).
        /// </summary>
        public static WaveletPack CDF26
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.353553390593274f,
                   0.707106781186548f,
                   0.353553390593274f,
                   }, new float[] {
                   0.006905339660025f,
                   0.013810679320050f,
                  -0.046956309688169f,
                  -0.107723298696388f,
                   0.169871355636612f,
                   0.447466009969612f,
                  -0.966747552403483f,
                   0.447466009969612f,
                   0.169871355636612f,
                  -0.107723298696388f,
                  -0.046956309688169f,
                   0.013810679320050f,
                   0.006905339660025f,
                   0.000000000000000f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 4/2).
        /// </summary>
        public static WaveletPack CDF42
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.088388347648318f,
                   0.353553390593274f,
                   0.530330085889911f,
                   0.353553390593274f,
                   0.088388347648318f
                                   }, new float[] {
                  -0.132582521472478f,
                  -0.530330085889911f,
                  -0.220970869120796f,
                   1.767766952966369f,
                  -0.220970869120796f,
                  -0.530330085889911f,
                  -0.132582521472478f,
                   0.000000000000000f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 4/4).
        /// </summary>
        public static WaveletPack CDF44
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.088388347648318f,
                   0.353553390593274f,
                   0.530330085889911f,
                   0.353553390593274f,
                   0.088388347648318f
                                   }, new float[] {
                   0.027621358640100f,
                   0.110485434560398f,
                   0.005524271728020f,
                  -0.530330085889911f,
                  -0.386699020961393f,
                   1.546796083845573f,
                  -0.386699020961393f,
                  -0.530330085889911f,
                   0.005524271728020f,
                   0.110485434560398f,
                   0.027621358640100f,
                   0.000000000000000f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 4/6).
        /// </summary>
        public static WaveletPack CDF46
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.088388347648318f,
                   0.353553390593274f,
                   0.530330085889911f,
                   0.353553390593274f,
                   0.088388347648318f
                                   }, new float[] {
                  -0.006042172202522f,
                  -0.024168688810087f,
                   0.009494842032534f,
                   0.158822812180572f,
                   0.096156854765846f,
                  -0.506161397079824f,
                  -0.453162915189133f,
                   1.450121328605225f,
                  -0.453162915189133f,
                  -0.506161397079824f,
                   0.096156854765846f,
                   0.158822812180572f,
                   0.009494842032534f,
                  -0.024168688810087f,
                  -0.006042172202522f,
                   0.000000000000000f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 6/2).
        /// </summary>
        public static WaveletPack CDF62
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.022097086912080f,
                   0.132582521472478f,
                   0.331456303681194f,
                   0.441941738241592f,
                   0.331456303681194f,
                   0.132582521472478f,
                   0.022097086912080f }, new float[] {
                   0.110485434560398f,
                   0.662912607362388f,
                   1.237436867076458f,
                  -0.309359216769115f,
                  -3.402951384460260f,
                  -0.309359216769115f,
                   1.237436867076458f,
                   0.662912607362388f,
                   0.110485434560398f,
                   0.000000000000000f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 6/4).
        /// </summary>
        public static WaveletPack CDF64
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.022097086912080f,
                   0.132582521472478f,
                   0.331456303681194f,
                   0.441941738241592f,
                   0.331456303681194f,
                   0.132582521472478f,
                   0.022097086912080f }, new float[] {
                  -0.024168688810087f,
                  -0.145012132860522f,
                  -0.227876208780821f,
                   0.324550964021169f,
                   1.261605555886545f,
                   0.174014559432627f,
                  -2.726228097777822f,
                   0.174014559432627f,
                   1.261605555886545f,
                   0.324550964021169f,
                  -0.227876208780821f,
                  -0.145012132860522f,
                  -0.024168688810087f,
                   0.000000000000000f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 6/6).
        /// </summary>
        public static WaveletPack CDF66
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                   0.022097086912080f,
                   0.132582521472478f,
                   0.331456303681194f,
                   0.441941738241592f,
                   0.331456303681194f,
                   0.132582521472478f,
                   0.022097086912080f }, new float[] {
                   0.005437954982270f,
                   0.032627729893618f,
                   0.041086770977148f,
                  -0.134136222895983f,
                  -0.380138948284370f,
                   0.096156854765846f,
                   1.196350096099310f,
                   0.358905028829793f,
                  -2.432578528735264f,
                   0.358905028829793f,
                   1.196350096099310f,
                   0.096156854765846f,
                  -0.380138948284370f,
                  -0.134136222895983f,
                   0.041086770977148f,
                   0.032627729893618f,
                   0.005437954982270f,
                   0.000000000000000f
                });
            }
        }
        /// <summary>
        /// Returns Cohen-Daubechies-Feaveau wavelet (CDF 9/7).
        /// </summary>
        public static WaveletPack CDF97
        {
            get
            {
                // Cohen–Daubechies–Feauveau wavelet:
                return WaveletPack.Create(new float[] {
                 3.782845550750114e-02f,
                -2.384946501955685e-02f,
                -1.106244044092826e-01f,
                 3.774028556128305e-01f,
                 8.526986789091245e-01f,
                 3.774028557909638e-01f,
                -1.106244045129673e-01f,
                -2.384946502829822e-02f,
                 3.782845552136610e-02f }, new float[] {
                 6.453888262876165e-02f,
                -4.068941760920477e-02f,
                -4.180922732220352e-01f,
                 7.884856164063713e-01f,
                -4.180922732220352e-01f,
                -4.068941760920475e-02f,
                 6.453888262876159e-02f,
                 0.000000000000000e-00f,
                 0.000000000000000e-00f
                });
            }
        }
        #endregion
    }
}
