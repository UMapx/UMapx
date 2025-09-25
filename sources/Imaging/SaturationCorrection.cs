using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Colorspace;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the saturation correction filter.
    /// </summary>
    [Serializable]
    public class SaturationCorrection : IBitmapFilter
    {
        #region Private data
        private float saturation;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the saturation correction filter.
        /// </summary>
        /// <param name="saturation">Saturation [-100, 100]</param>
        public SaturationCorrection(float saturation)
        {
            Saturation = saturation;
        }
        /// <summary>
        /// Initializes the saturation correction filter.
        /// </summary>
        public SaturationCorrection()
        {
            Saturation = 20;
        }
        /// <summary>
        /// Gets or sets the saturation value [-100, 100].
        /// </summary>
        public float Saturation
        {
            get
            {
                return this.saturation;
            }
            set
            {
                this.saturation = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            if (bmData.PixelFormat != PixelFormat.Format32bppArgb)
                throw new NotSupportedException("Only support Format32bppArgb pixelFormat");

            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            float s = this.saturation / 100.0f;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;

                    byte r = p[k + 2];
                    byte g = p[k + 1];
                    byte b = p[k + 0];

                    RGB rgb = RGB.Saturation(r, g, b, s);

                    p[k + 2] = rgb.Red;
                    p[k + 1] = rgb.Green;
                    p[k + 0] = rgb.Blue;
                }
            });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            Apply(bmData);
            BitmapFormat.Unlock(Data, bmData);
        }
        #endregion
    }
}
