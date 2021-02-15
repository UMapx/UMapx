using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the local contrast inversion filter.
    /// <remarks>
    /// This filter is used to equalize the illumination of images by averaging the brightness.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LocalContrastInversion : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private float a;
        private float b;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local contrast inversion filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public LocalContrastInversion(int radius, Space space, float a = 0.75f, float b = 0.05f)
        {
            this.gb = new BoxBlur(radius);
            Space = space;
            A = a;
            B = b;
        }
        /// <summary>
        /// Initializes the local contrast inversion filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public LocalContrastInversion(int width, int height, Space space, float a = 0.75f, float b = 0.05f)
        {
            this.gb = new BoxBlur(width, height);
            Space = space;
            A = a;
            B = b;
        }
        /// <summary>
        /// Initializes the local contrast inversion filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public LocalContrastInversion(SizeInt size, Space space, float a = 0.75f, float b = 0.05f)
        {
            this.gb = new BoxBlur(size);
            Space = space;
            A = a;
            B = b;
        }
        /// <summary>
        /// Gets or sets the contrast value (0, 1].
        /// </summary>
        public float A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the offset value (0, 1].
        /// </summary>
        public float B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.b = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.LocalContrastInversion(this.a, this.b, 256);
        }
        #endregion
    }
}
