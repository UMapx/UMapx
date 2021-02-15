using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the filter for homomorphic processing.
    /// <remarks>
    /// A homomorphic filter is most often used to equalize the illumination of images.
    /// It simultaneously normalizes the brightness of the image and increases the contrast.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Homomorphic_filtering
    /// </remarks>
    /// </summary>
    [Serializable]
    public class HomomorphicEnhancement : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private float a;
        private float b;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the filter for homomorphic processing.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public HomomorphicEnhancement(int radius, Space space, float a = 0.5f, float b = 0.05f)
        {
            gb = new BoxBlur(radius);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Initializes the filter for homomorphic processing.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public HomomorphicEnhancement(int width, int height, Space space, float a = 0.5f, float b = 0.05f)
        {
            gb = new BoxBlur(width, height);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Initializes the filter for homomorphic processing.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public HomomorphicEnhancement(SizeInt size, Space space, float a = 0.5f, float b = 0.05f)
        {
            gb = new BoxBlur(size);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the contrast value [-1, 1].
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
            this.values = Intensity.HomomorphicEnhancement(this.a, this.b, 256);
        }
        #endregion
    }
}
