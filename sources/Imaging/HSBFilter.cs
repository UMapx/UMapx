using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the HSB filter.
    /// </summary>
    [Serializable]
    public class HSBFilter : IBitmapFilter
    {
        #region Private data
        private int hue;
        private float saturation;
        private float brightness;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the HSB filter.
        /// </summary>
        /// <param name="hue">Hue [-180, 180]</param>
        /// <param name="saturation">Saturation [-1, 1]</param>
        /// <param name="brightness">Brightness [-1, 1]</param>
        public HSBFilter(int hue, float saturation, float brightness)
        {
            Hue = hue;
            Saturation = saturation;
            Brightness = brightness;
        }
        /// <summary>
        /// Initializes the HSB filter.
        /// </summary>
        public HSBFilter()
        {
            new HSBFilter(0, 0, 0);
        }
        /// <summary>
        /// Hue [-180, 180].
        /// </summary>
        public int Hue
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
        /// Brightness [-1, 1].
        /// </summary>
        public float Brightness
        {
            get
            {
                return this.brightness;
            }
            set
            {
                this.brightness = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                HSB hsb; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    hsb = HSB.FromRGB(p[k + 2], p[k + 1], p[k + 0]);

                    hsb.Hue = Maths.Scale(hsb.Hue + hue, 0, 360);
                    hsb.Saturation += saturation;
                    hsb.Brightness += brightness;
                    rgb = hsb.ToRGB;

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
            return;
        }
        #endregion
    }
}
