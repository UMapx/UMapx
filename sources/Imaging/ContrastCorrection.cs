using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the contrast correction filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://esate.ru/uroki/OpenGL/image_processing/_p4106/
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ContrastCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private double contrast;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the contrast correction filter.
        /// </summary>
        /// <param name="value">Contrast [-1, 1]</param>
        /// <param name="space">Color space</param>
        public ContrastCorrection(double value, Space space)
        {
            Contrast = value; this.Space = space;
        }
        /// <summary>
        /// Initializes the contrast correction filter.
        /// </summary>
        public ContrastCorrection()
        {
            Contrast = 0.5;
        }
        /// <summary>
        /// Gets or sets the contrast value.
        /// </summary>
        public double Contrast
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
            this.values = Intensity.Contrast(this.contrast, 256);
        }
        #endregion
    }
}
