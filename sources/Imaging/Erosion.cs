using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the erosion filter.
    /// </summary>
    [Serializable]
    public class Erosion : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public Erosion(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public Erosion(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public Erosion(SizeInt size)
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
                ApplyHorizontal(bmSrc, bmData);
                ApplyVertical(bmData, bmSrc);
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
            int r = rh >> 1;
            #endregion

            Parallel.For(0, height, y =>
            {
                byte red, green, blue;
                int x, i, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y - r;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x;
                    v = ystride + x * 4;

                    red = green = blue = 255;

                    #region Convolution filtering
                    for (i = 0; i < rh; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        jr = xr + 0;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red = p[2] < red ? p[2] : red;
                        green = p[1] < green ? p[1] : green;
                        blue = p[0] < blue ? p[0] : blue;
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red;
                    dst[v + 1] = green;
                    dst[v] = blue;
                    #endregion
                }
            }
            );

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
            int r = rw >> 1;
            #endregion

            Parallel.For(0, height, y =>
            {
                byte red, green, blue;
                int x, j, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x - r;
                    v = ystride + x * 4;

                    red = green = blue = 255;

                    #region Convolution filtering
                    ir = yr + 0;
                    if (ir < 0) continue; if (ir >= height) break;

                    for (j = 0; j < rw; j++)
                    {
                        jr = xr + j;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red = p[2] < red ? p[2] : red;
                        green = p[1] < green ? p[1] : green;
                        blue = p[0] < blue ? p[0] : blue;
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red;
                    dst[v + 1] = green;
                    dst[v] = blue;
                    #endregion
                }
            }
            );

            return;
        }
        #endregion
    }
}
