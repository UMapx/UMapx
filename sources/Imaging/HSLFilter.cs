using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the HSL filter.
    /// </summary>
    [Serializable]
    public class HSLFilter : IBitmapFilter
    {
        #region Private data
        private float hue;
        private float saturation;
        private float lightness;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the HSL filter.
        /// </summary>
        /// <param name="hue">Hue [-180, 180]</param>
        /// <param name="saturation">Saturation [-1, 1]</param>
        /// <param name="lightness">Lightness [-1, 1]</param>
        public HSLFilter(float hue, float saturation, float lightness)
        {
            Hue = hue;
            Saturation = saturation;
            Lightness = lightness;
        }
        /// <summary>
        /// Initializes the HSL filter.
        /// </summary>
        public HSLFilter()
        {
            new HSLFilter(0, 0, 0);
        }
        /// <summary>
        /// Hue [-180, 180].
        /// </summary>
        public float Hue
        {
            get
            {
                return this.hue;
            }
            set
            {
                this.hue = Maths.Scale(value, 0, 360);
            }
        }
        /// <summary>
        /// Saturation [-1, 1].
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
        /// Lightness [-1, 1].
        /// </summary>
        public float Lightness
        {
            get
            {
                return this.lightness;
            }
            set
            {
                this.lightness = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                HSL hsl; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    hsl = HSL.FromRGB(p[k + 2], p[k + 1], p[k + 0]);

                    hsl.Hue = Maths.Scale(hsl.Hue + hue, 0, 360);
                    hsl.Saturation += saturation;
                    hsl.Lightness += lightness;
                    rgb = hsl.ToRGB;

                    p[k + 0] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
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
