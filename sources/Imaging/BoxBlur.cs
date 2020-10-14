using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the box blur filter.
    /// </summary>
    [Serializable]
    public class BoxBlur : IBitmapFilter
    {
        #region Private data
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the box blur filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public BoxBlur(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Initializes the box blur filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public BoxBlur(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the box blur filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public BoxBlur(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(rw, rh);
            }
            set
            {
                this.rw = value.Width;
                this.rh = value.Height;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            if (rw >= 2 && rh >= 2)
            {
                ApplyHorizontal(bmData);
                ApplyVertical(bmData);
            }
            else if (rw >= 2 && rh < 2)
            {
                ApplyHorizontal(bmData);
            }
            else if (rw < 2 && rh >= 2)
            {
                ApplyVertical(bmData);
            }
            else return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyVertical(BitmapData bmData)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int h = rh >= height ? height - 1 : rh;
            int v = h >> 1;
            int dl = height - v;
            #endregion

            Parallel.For(0, width, x =>
            {
                double r = 0;
                double g = 0;
                double b = 0;
                int p, w, q, y;
                int xx = x * 4;

                for (p = xx, y = 0; y < h; y++, p += stride)
                {
                    r += dst[p + 2];
                    g += dst[p + 1];
                    b += dst[p + 0];
                }

                for (p = xx, y = 0; y < v; y++, p += stride)
                {
                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

                for (
                    y = v,
                    p = xx,
                    q = xx + (y + 0) * stride,
                    w = xx + (y + v) * stride;

                    y < dl;

                    y++,
                    p += stride,
                    q += stride,
                    w += stride)
                {
                    r = r - dst[p + 2] + dst[w + 2];
                    g = g - dst[p + 1] + dst[w + 1];
                    b = b - dst[p + 0] + dst[w + 0];

                    dst[q + 2] = Maths.Byte(r / h);
                    dst[q + 1] = Maths.Byte(g / h);
                    dst[q + 0] = Maths.Byte(b / h);
                }

                for (
                    y = dl,
                    w = xx + (y - v) * stride,
                    p = xx + (y + 0) * stride;

                    y < height;

                    y++,
                    w += stride,
                    p += stride)
                {
                    r = r - dst[w + 2] + dst[p + 2];
                    g = g - dst[w + 1] + dst[p + 1];
                    b = b - dst[w + 0] + dst[p + 0];

                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

            });

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyHorizontal(BitmapData bmData)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int h = rw >= width ? width - 1 : rw;
            int v = h >> 1;
            int dl = width - v;
            #endregion

            Parallel.For(0, height, y =>
            {
                double r = 0;
                double g = 0;
                double b = 0;
                int p, q, w, x;
                int yy = y * stride;

                for (p = yy, x = 0; x < h; x++, p += 4)
                {
                    r += dst[p + 2];
                    g += dst[p + 1];
                    b += dst[p + 0];
                }

                for (p = yy, x = 0; x < v; x++, p += 4)
                {
                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

                for (
                    x = v,
                    p = yy,
                    q = yy + (x + 0) * 4,
                    w = yy + (x + v) * 4;

                    x < dl;

                    x++,
                    p += 4,
                    q += 4,
                    w += 4)
                {
                    r = r - dst[p + 2] + dst[w + 2];
                    g = g - dst[p + 1] + dst[w + 1];
                    b = b - dst[p + 0] + dst[w + 0];

                    dst[q + 2] = Maths.Byte(r / h);
                    dst[q + 1] = Maths.Byte(g / h);
                    dst[q + 0] = Maths.Byte(b / h);
                }

                for (
                    x = dl,
                    w = (x - v) * 4 + yy,
                    p = (x + 0) * 4 + yy;

                    x < width;

                    x++,
                    p += 4,
                    w += 4)
                {
                    r = r - dst[w + 2] + dst[p + 2];
                    g = g - dst[w + 1] + dst[p + 1];
                    b = b - dst[w + 0] + dst[p + 0];

                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

            });
        }
        #endregion
    }
}
