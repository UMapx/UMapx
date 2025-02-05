using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the box blur filter.
    /// </summary>
    [Serializable]
    public class BoxBlur : IBitmapFilter2, IBitmapFilter
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
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            if (rw >= 2 && rh >= 2)
            {
                ApplyVertical(bmSrc, bmData);
                ApplyHorizontal(bmData, bmSrc);
            }
            else if (rw >= 2 && rh < 2)
            {
                ApplyHorizontal(bmData, bmSrc);
            }
            else if (rw < 2 && rh >= 2)
            {
                ApplyVertical(bmData, bmSrc);
            }
            else return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Data, bmData);
            BitmapFormat.Unlock(Src, bmSrc);
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            Bitmap current = BitmapFormat.Bitmap(bmData);
            Bitmap Src = (Bitmap)current.Clone();
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Src, bmSrc);
            Src.Dispose();
            current.Dispose();
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            var Src = (Bitmap)Data.Clone();
            Apply(Data, Src);
            Src.Dispose();
            return;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyVertical(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int h = rh >= height ? height - 1 : rh;
            int v = h >> 1;
            int dl = height - v;
            #endregion

            Parallel.For(0, width, x =>
            {
                float r = 0;
                float g = 0;
                float b = 0;
                int p, w, q, y;
                int xx = x * 4;

                for (p = xx, y = 0; y < h; y++, p += stride)
                {
                    r += src[p + 2];
                    g += src[p + 1];
                    b += src[p + 0];
                }

                for (p = xx, y = 0; y < v; y++, p += stride)
                {
                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

                for (
                    y = v,
                    p = xx + (y - v) * stride,
                    q = xx + (y + 0) * stride,
                    w = xx + (y + v) * stride;

                    y < dl;

                    y++,
                    p += stride,
                    q += stride,
                    w += stride)
                {
                    r = r - src[p + 2] + src[w + 2];
                    g = g - src[p + 1] + src[w + 1];
                    b = b - src[p + 0] + dst[w + 0];

                    dst[q + 2] = Maths.Byte(r / h);
                    dst[q + 1] = Maths.Byte(g / h);
                    dst[q + 0] = Maths.Byte(b / h);
                }

                for (
                    y = dl,
                    p = xx + (y - v) * stride,
                    q = xx + (y + 0) * stride;

                    y < height;

                    y++,
                    p += stride,
                    q += stride)
                {
                    r = r - src[p + 2] + src[q + 2];
                    g = g - src[p + 1] + src[q + 1];
                    b = b - src[p + 0] + src[q + 0];

                    dst[q + 2] = Maths.Byte(r / h);
                    dst[q + 1] = Maths.Byte(g / h);
                    dst[q + 0] = Maths.Byte(b / h);
                }

            });

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyHorizontal(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int h = rw >= width ? width - 1 : rw;
            int v = h >> 1;
            int dl = width - v;
            #endregion

            Parallel.For(0, height, y =>
            {
                float r = 0;
                float g = 0;
                float b = 0;
                int p, q, w, x;
                int yy = y * stride;

                for (p = yy, x = 0; x < h; x++, p += 4)
                {
                    r += src[p + 2];
                    g += src[p + 1];
                    b += src[p + 0];
                }

                for (p = yy, x = 0; x < v; x++, p += 4)
                {
                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

                for (
                    x = v,
                    p = yy + (x - v) * 4,
                    q = yy + (x + 0) * 4,
                    w = yy + (x + v) * 4;

                    x < dl;

                    x++,
                    p += 4,
                    q += 4,
                    w += 4)
                {
                    r = r - src[p + 2] + src[w + 2];
                    g = g - src[p + 1] + src[w + 1];
                    b = b - src[p + 0] + src[w + 0];

                    dst[q + 2] = Maths.Byte(r / h);
                    dst[q + 1] = Maths.Byte(g / h);
                    dst[q + 0] = Maths.Byte(b / h);
                }

                for (
                    x = dl,
                    p = (x - v) * 4 + yy,
                    q = (x + 0) * 4 + yy;

                    x < width;

                    x++,
                    p += 4,
                    q += 4)
                {
                    r = r - src[p + 2] + src[q + 2];
                    g = g - src[p + 1] + src[q + 1];
                    b = b - src[p + 0] + src[q + 0];

                    dst[q + 2] = Maths.Byte(r / h);
                    dst[q + 1] = Maths.Byte(g / h);
                    dst[q + 0] = Maths.Byte(b / h);
                }

            });
        }
        #endregion
    }
}
