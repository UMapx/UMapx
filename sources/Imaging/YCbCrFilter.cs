using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the YCbCr filter.
    /// </summary>
    [Serializable]
    public class YCbCrFilter : IBitmapFilter
    {
        #region Private data
        private float y;
        private float cb;
        private float cr;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the YCbCr filter.
        /// </summary>
        /// <param name="y">Y [-1, 1]</param>
        /// <param name="cb">Cb [-1, 1]</param>
        /// <param name="cr">Cr [-1, 1]</param>
        public YCbCrFilter(float y, float cb, float cr)
        {
            Y = y;
            Cb = cb;
            Cr = cr;
        }
        /// <summary>
        /// Initializes the YCbCr filter.
        /// </summary>
        public YCbCrFilter()
        {
            new YCbCrFilter(0, 0, 0);
        }
        /// <summary>
        /// Y [-1, 1].
        /// </summary>
        public float Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = value;
            }
        }
        /// <summary>
        /// Cb [-1, 1].
        /// </summary>
        public float Cb
        {
            get
            {
                return this.cb;
            }
            set
            {
                this.cb = value;
            }
        }
        /// <summary>
        /// Cr [-1, 1].
        /// </summary>
        public float Cr
        {
            get
            {
                return this.cr;
            }
            set
            {
                this.cr = value;
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
                YCbCr ycbcr; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    ycbcr = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    ycbcr.Y += y;
                    ycbcr.Cb += cb;
                    ycbcr.Cr += cr;
                    rgb = ycbcr.ToRGB;

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
