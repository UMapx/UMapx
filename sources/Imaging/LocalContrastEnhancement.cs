using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the local contrast enhancement filter.
    /// <remarks>
    /// This filter is also known as "Unsharp Masking."
    /// More information can be found on the website:
    /// http://www.cambridgeincolour.com/tutorials/local-contrast-enhancement.htm
    /// Filter usage example:
    /// http://www.knowhowtransfer.com/photoshop-professional-plugins/alce-local-contrast-enhancer/
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LocalContrastEnhancement : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private float contrast;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local contrast enhancement filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [-1, 1]</param>
        public LocalContrastEnhancement(int radius, Space space, float contrast = 0.75f)
        {
            gb = new BoxBlur(radius);
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Initializes the local contrast enhancement filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [-1, 1]</param>
        public LocalContrastEnhancement(int width, int height, Space space, float contrast = 0.75f)
        {
            gb = new BoxBlur(width, height);
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Initializes the local contrast enhancement filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [-1, 1]</param>
        public LocalContrastEnhancement(SizeInt size, Space space, float contrast = 0.75f)
        {
            gb = new BoxBlur(size);
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Gets or sets the contrast value [-1, 1].
        /// </summary>
        public float Contrast
        {
            get
            {
                return this.contrast;
            }
            set
            {
                this.contrast = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.LocalContrastEnhancement(this.contrast, 256);
        }
        #endregion
    }
}
