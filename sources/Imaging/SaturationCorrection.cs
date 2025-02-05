using System;
using SkiaDrawing;
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
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            float s = this.saturation / 255.0f;

            Parallel.For(0, height, y =>
            {
                RGB rgb;
                int x, ystride, k;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    rgb = RGB.Saturation(p[k + 2], p[k + 1], p[k], s);

                    p[k] = rgb.Red;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Blue;
                }
            }
            );

            return;
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
