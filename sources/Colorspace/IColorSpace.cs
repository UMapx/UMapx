namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines the color space interface.
    /// </summary>
    public interface IColorSpace
    {
        #region Interface
        /// <summary>
        /// Returns the color model RGB.
        /// </summary>
        RGB ToRGB { get; }
        #endregion
    }
}
