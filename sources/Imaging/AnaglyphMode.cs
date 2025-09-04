namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the stereo effect creation algorithm.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.3dtv.at/Knowhow/AnaglyphComparison_en.aspx
    /// </remarks>
    public enum AnaglyphMode
    {
        /// <summary>
        /// Create a stereo effect for a pair of images using the following calculations:
        /// <list type="bullet">
        /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
        /// <item>G<sub>a</sub>=0;</item>
        /// <item>B<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>.</item>
        /// </list>
        /// </summary>
        True,
        /// <summary>
        /// Create a stereo effect for a pair of images using the following calculations:
        /// <list type="bullet">
        /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
        /// <item>G<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>;</item>
        /// <item>B<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>.</item>
        /// </list>
        /// </summary>
        Gray,
        /// <summary>
        /// Create a stereo effect for a pair of images using the following calculations:
        /// <list type="bullet">
        /// <item>R<sub>a</sub>=R<sub>l</sub>;</item>
        /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
        /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
        /// </list>
        /// </summary>
        Color,
        /// <summary>
        /// Create a stereo effect for a pair of images using the following calculations:
        /// <list type="bullet">
        /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
        /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
        /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
        /// </list>
        /// </summary>
        HalfColor,
        /// <summary>
        /// Create a stereo effect for a pair of images using the following calculations:
        /// <list type="bullet">
        /// <item>R<sub>a</sub>=0.7*G<sub>l</sub>+0.3*B<sub>l</sub>;</item>
        /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
        /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
        /// </list>
        /// </summary>
        Optimized
    }
}
