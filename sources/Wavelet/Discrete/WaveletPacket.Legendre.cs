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
        #region Legendre wavelets
        /// <summary>
        /// Returns Legendre wavelet of 1 order.
        /// </summary>
        /// <remarks>
        /// Haar wavelet.
        /// </remarks>
        public static WaveletPacket L1
        {
            get
            {
                // Haar's wavelet:
                return WaveletPacket.Haar;
            }
        }
        /// <summary>
        /// Returns Legendre wavelet of 2 order.
        /// </summary>
        /// <remarks>
        /// Nonorthogonal wavelet.
        /// </remarks>
        public static WaveletPacket L2
        {
            get
            {
                return Create(new float[] { 0.441941738001176f, 0.265165042944955f, 0.265165042944955f, 0.441941738001176f });
            }
        }
        /// <summary>
        /// Returns Legendre wavelet of 3 order.
        /// </summary>
        /// <remarks>
        /// Nonorthogonal wavelet.
        /// </remarks>
        public static WaveletPacket L3
        {
            get
            {
                return Create(new float[] { 0.348029118865254f, 0.193349510480697f, 0.165728151840597f, 0.165728151840597f, 0.193349510480697f, 0.348029118865254f });
            }
        }
        /// <summary>
        /// Returns Legendre wavelet of 4 order.
        /// </summary>
        /// <remarks>
        /// Nonorthogonal wavelet.
        /// </remarks>
        public static WaveletPacket L4
        {
            get
            {
                return Create(new float[] {
                0.209472656663610f,
                0.112792968646356f,
                0.092285156550895f,
                0.085449218714337f,
                0.085449218714337f,
                0.092285156550895f,
                0.112792968646356f,
                0.209472656663610f });
            }
        }
        /// <summary>
        /// Returns Legendre wavelet of 5 order.
        /// </summary>
        /// <remarks>
        /// Nonorthogonal wavelet.
        /// </remarks>
        public static WaveletPacket L5
        {
            get
            {
                return Create(new float[] {
                0.185470581656607f,
                0.098190308518174f,
                0.078552246248854f,
                0.070495604520089f,
                0.067291259631474f,
                0.067291259631474f,
                0.070495604520089f,
                0.078552246248854f,
                0.098190308518174f,
                0.185470581656607f });
            }
        }
        /// <summary>
        /// Returns Legendre wavelet of 6 order.
        /// </summary>
        /// <remarks>
        /// Nonorthogonal wavelet.
        /// </remarks>
        public static WaveletPacket L6
        {
            get
            {
                return Create(new float[] {
                0.204303513760092f,
                0.105674265879069f,
                0.082191094112372f,
                0.071232282506865f,
                0.065038168525027f,
                0.061321700135924f,
                0.059170059053587f,
                0.058175612360798f,
                0.058175612360798f,
                0.059170059053587f,
                0.061321700135924f,
                0.065038168525027f,
                0.071232282506865f,
                0.082191094112372f,
                0.105674265879069f,
                0.204303513760092f });
            }
        }
        /// <summary>
        /// Returns Legendre wavelet of 7 order.
        /// </summary>
        /// <remarks>
        /// Nonorthogonal wavelet.
        /// </remarks>
        public static WaveletPacket L7
        {
            get
            {
                return Create(new float[] {
                0.204303513760092f,
                0.105674265879069f,
                0.082191094112372f,
                0.071232282506865f,
                0.065038168525027f,
                0.061321700135924f,
                0.059170059053587f,
                0.058175612360798f,
                0.058175612360798f,
                0.059170059053587f,
                0.061321700135924f,
                0.065038168525027f,
                0.071232282506865f,
                0.082191094112372f,
                0.105674265879069f,
                0.204303513760092f
                });
            }
        }
        /// <summary>
        /// Returns Legendre wavelet of 8 order.
        /// </summary>
        /// <remarks>
        /// Nonorthogonal wavelet.
        /// </remarks>
        public static WaveletPacket L8
        {
            get
            {
                return Create(new float[] {
                0.192098002188675f,
                0.098959551600650f,
                0.076613845944306f,
                0.066046417942185f,
                0.059931004945093f,
                0.056095417347632f,
                0.053656492922234f,
                0.052196444698305f,
                0.051509660166010f,
                0.051509660166010f,
                0.052196444698305f,
                0.053656492922234f,
                0.056095417347632f,
                0.059931004945093f,
                0.066046417942185f,
                0.076613845944306f,
                0.098959551600650f,
                0.192098002188675f
                });
            }
        }
        /// <summary>
        /// Returns Legendre wavelet of 9 order.
        /// </summary>
        /// <remarks>
        /// Nonorthogonal wavelet.
        /// </remarks>
        public static WaveletPacket L9
        {
            get
            {
                return Create(new float[] {
                0.181847075181813f,
                0.093380945787564f,
                0.072036729607550f,
                0.061849710911572f,
                0.055992816153682f,
                0.052011550417160f,
                0.049443084029449f,
                0.047747894516504f,
                0.046709890045993f,
                0.046215608263808f,
                0.046215608263808f,
                0.046709890045993f,
                0.047747894516504f,
                0.049443084029449f,
                0.052011550417160f,
                0.055992816153682f,
                0.061849710911572f,
                0.072036729607550f,
                0.093380945787564f,
                0.181847075181813f
                });
            }
        }
        #endregion
    }
}
