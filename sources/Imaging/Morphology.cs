using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the morphology filter.
    /// </summary>
    [Serializable]
    public class Morphology : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int rw;
        private int rh;
        private int tw;
        private int th;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="threshold">Threshold</param>
        public Morphology(int radius = 3, int threshold = 0)
        {
            Size = new SizeInt(radius, radius);
            Threshold = new SizeInt(threshold, threshold);
        }
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="widthThreshold">Threshold by width</param>
        /// <param name="heightThreshold">Threshold by height</param>
        public Morphology(int width, int height, int widthThreshold, int heightThreshold)
        {
            Size = new SizeInt(width, height);
            Threshold = new SizeInt(widthThreshold, heightThreshold);
        }
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="threshold">Thresholds</param>
        public Morphology(SizeInt size, SizeInt threshold)
        {
            Size = size;
            Threshold = threshold;
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
        /// Gets or sets thresholds.
        /// </summary>
        public SizeInt Threshold
        {
            get
            {
                return new SizeInt(tw, th);
            }
            set
            {
                this.tw = value.Width;
                this.th = value.Height;
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
            int threshold = th;
            #endregion

            Parallel.For(0, width, x =>
            {
                byte[] r = new byte[h];
                byte[] g = new byte[h];
                byte[] b = new byte[h];
                int p, w, q, y, i;
                int xx = x * 4;

                for (p = xx, y = 0; y < h; y++, p += stride)
                {
                    r[y] = src[p + 2];
                    g[y] = src[p + 1];
                    b[y] = src[p + 0];
                }

                Array.Sort(r);
                Array.Sort(g);
                Array.Sort(b);

                for (p = xx, y = 0; y < v; y++, p += stride)
                {
                    dst[p + 2] = r[threshold];
                    dst[p + 1] = g[threshold];
                    dst[p + 0] = b[threshold];
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
                    i = Array.IndexOf(r, src[p + 2]);
                    r[i] = src[w + 2];
                    FastSort(ref r, i);
                    dst[q + 2] = r[threshold];

                    i = Array.IndexOf(g, src[p + 1]);
                    g[i] = src[w + 1];
                    FastSort(ref g, i);
                    dst[q + 1] = g[threshold];

                    i = Array.IndexOf(b, src[p + 0]);
                    b[i] = src[w + 0];
                    FastSort(ref b, i);
                    dst[q + 0] = b[threshold];
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
                    i = Array.IndexOf(r, src[p + 2]);
                    r[i] = src[q + 2];
                    FastSort(ref r, i);
                    dst[q + 2] = r[threshold];

                    i = Array.IndexOf(g, src[p + 1]);
                    g[i] = src[q + 1];
                    FastSort(ref g, i);
                    dst[q + 1] = g[threshold];

                    i = Array.IndexOf(b, src[p + 0]);
                    b[i] = src[q + 0];
                    FastSort(ref b, i);
                    dst[q + 0] = b[threshold];
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
            int threshold = tw;
            #endregion

            Parallel.For(0, height, y =>
            {
                byte[] r = new byte[h];
                byte[] g = new byte[h];
                byte[] b = new byte[h];
                int p, q, w, x, i;
                int yy = y * stride;

                for (p = yy, x = 0; x < h; x++, p += 4)
                {
                    r[x] = src[p + 2];
                    g[x] = src[p + 1];
                    b[x] = src[p + 0];
                }

                Array.Sort(r);
                Array.Sort(g);
                Array.Sort(b);

                for (p = yy, x = 0; x < v; x++, p += 4)
                {
                    dst[p + 2] = r[threshold];
                    dst[p + 1] = g[threshold];
                    dst[p + 0] = b[threshold];
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
                    i = Array.IndexOf(r, src[p + 2]);
                    r[i] = src[w + 2];
                    FastSort(ref r, i);
                    dst[q + 2] = r[threshold];

                    i = Array.IndexOf(g, src[p + 1]);
                    g[i] = src[w + 1];
                    FastSort(ref g, i);
                    dst[q + 1] = g[threshold];

                    i = Array.IndexOf(b, src[p + 0]);
                    b[i] = src[w + 0];
                    FastSort(ref b, i);
                    dst[q + 0] = b[threshold];
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
                    i = Array.IndexOf(r, src[p + 2]);
                    r[i] = src[q + 2];
                    FastSort(ref r, i);
                    dst[q + 2] = r[threshold];

                    i = Array.IndexOf(g, src[p + 1]);
                    g[i] = src[q + 1];
                    FastSort(ref g, i);
                    dst[q + 1] = g[threshold];

                    i = Array.IndexOf(b, src[p + 0]);
                    b[i] = src[q + 0];
                    FastSort(ref b, i);
                    dst[q + 0] = b[threshold];
                }

            });
        }
        /// <summary>
        /// O(N) sort algorithm.
        /// </summary>
        /// <param name="s">Array</param>
        /// <param name="index">Index</param>
        private static void FastSort(ref byte[] s, int index)
        {
            int length = s.Length - 1;

            for (int i = index; i < length; i++)
            {
                if (s[i] > s[i + 1])
                {
                    var t = s[i + 1];
                    s[i + 1] = s[i];
                    s[i] = t;
                }
                else
                    break;
            }

            for (int i = index; i > 0; i--)
            {
                if (s[i] < s[i - 1])
                {
                    var t = s[i - 1];
                    s[i - 1] = s[i];
                    s[i] = t;
                }
                else
                    break;
            }
        }
        #endregion

        #region Public static
        /// <summary>
        /// Returns erosion filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <returns>Morphology</returns>
        public static Morphology Erosion(int width, int height)
        {
            return new Morphology(width, height, 0, 0);
        }
        /// <summary>
        /// Returns dilatation filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <returns>Morphology</returns>
        public static Morphology Dilatation(int width, int height)
        {
            return new Morphology(width, height, width - 1, height - 1);
        }
        /// <summary>
        /// Returns median filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <returns>Morphology</returns>
        public static Morphology Median(int width, int height)
        {
            return new Morphology(width, height, width / 2, height / 2);
        }
        #endregion
    }
}
