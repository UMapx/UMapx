using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the global contrast enhancement filter.
    /// </summary>
    [Serializable]
    public class ContrastEnhancement : Correction, IBitmapFilter
    {
        #region Private data
        private float contrast;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the global contrast enhancement filter.
        /// </summary>
        /// <param name="contrast">Contrast [-1, 1]</param>
        /// <param name="space">Color space</param>
        public ContrastEnhancement(float contrast, Space space)
        {
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Gets or sets the contrast coefficent value [-1, 1].
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
            this.values = Intensity.LogContrast(1 + this.contrast, 256);
        }
        #endregion
    }
}
