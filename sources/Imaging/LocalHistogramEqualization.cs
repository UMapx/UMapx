using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the local histogram equalization filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://angeljohnsy.blogspot.com/2011/06/local-histogram-equalization.html
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LocalHistogramEqualization : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int l0;
        private int l1;
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local histogram equalization filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public LocalHistogramEqualization(int radius = 10)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Initializes the local histogram equalization filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public LocalHistogramEqualization(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the local histogram equalization filter.
        /// </summary>
        /// <param name="size">Radius</param>
        public LocalHistogramEqualization(SizeInt size)
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
                return new SizeInt(l0, l1);
            }
            set
            {
                Data(value);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="size">Radius</param>
        private void Data(SizeInt size)
        {
            this.l0 = size.Width;
            this.l1 = size.Height;
            this.rw = l0 >> 1;
            this.rh = l1 >> 1;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            #endregion

            Parallel.For(0, height, y =>
            {
                int[] H = new int[256];
                int[] cdf = new int[256];
                int n = l0 * l1;
                int brightness;
                float dn = 255.0f / n;

                int x, i, j, ir, jr, yr, xr, irstride;
                int ystride, v;
                byte* p;

                yr = y - rw;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x - rh;
                    v = ystride + x * 4;

                    #region Convolution filtering
                    for (i = 0; i < l0; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;
                        irstride = ir * stride;

                        for (j = 0; j < l1; j++)
                        {
                            jr = xr + j;

                            if (jr < 0) continue; if (jr >= width) break;

                            p = &src[irstride + jr * 4];
                            brightness = (p[2] + p[1] + p[0]) / 3;
                            H[brightness]++;
                        }
                    }
                    #endregion

                    #region Density function
                    cdf[0] = H[0];

                    for (i = 1; i < 256; i++)
                    {
                        cdf[i] = H[i] + cdf[i - 1];
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = Maths.Byte(cdf[src[v + 2]] * dn);
                    dst[v + 1] = Maths.Byte(cdf[src[v + 1]] * dn);
                    dst[v] = Maths.Byte(cdf[src[v]] * dn);
                    #endregion

                    #region Clear data
                    Array.Clear(H, 0, 256);
                    Array.Clear(cdf, 0, 256);
                    #endregion
                }
            }
            );

            return;
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
    }
}
