using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the morphology filter.
    /// </summary>
    [Serializable]
    public class Morphology : IBitmapFilter2
    {
        #region Private data
        private int threshold;
        private int rw;
        private int rh;
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
            Threshold = threshold;
        }
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="threshold">Threshold</param>
        public Morphology(int width, int height, int threshold)
        {
            Size = new SizeInt(width, height);
            Threshold = threshold;
        }
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="threshold">Threshold</param>
        public Morphology(SizeInt size, int threshold = 0)
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
        /// Gets or sets the threshold.
        /// </summary>
        public int Threshold
        {
            get
            {
                return this.threshold;
            }
            set
            {
                this.threshold = value;
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
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Data, bmData);
            BitmapFormat.Unlock(Src, bmSrc);
            Src.Dispose();
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
                byte[] red, green, blue;
                int x, i, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y - r;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x;
                    v = ystride + x * 4;

                    red = new byte[rh];
                    green = new byte[rh];
                    blue = new byte[rh];

                    #region Convolution filtering
                    for (i = 0; i < rh; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        jr = xr + 0;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red[i] = p[2];
                        green[i] = p[1];
                        blue[i] = p[0];
                    }
                    #endregion

                    #region Morphology filtering
                    Array.Sort(red);
                    Array.Sort(green);
                    Array.Sort(blue);
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red[threshold];
                    dst[v + 1] = green[threshold];
                    dst[v] = blue[threshold];
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
                byte[] red, green, blue;
                int x, j, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x - r;
                    v = ystride + x * 4;

                    red = new byte[rw];
                    green = new byte[rw];
                    blue = new byte[rw];

                    #region Convolution filtering
                    ir = yr + 0;
                    if (ir < 0) continue; if (ir >= height) break;

                    for (j = 0; j < rw; j++)
                    {
                        jr = xr + j;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red[j] = p[2];
                        green[j] = p[1];
                        blue[j] = p[0];
                    }
                    #endregion

                    #region Morphology filtering
                    Array.Sort(red);
                    Array.Sort(green);
                    Array.Sort(blue);
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red[threshold];
                    dst[v + 1] = green[threshold];
                    dst[v] = blue[threshold];
                    #endregion
                }
            }
            );

            return;
        }
        #endregion

        #region Public static components
        /// <summary>
        /// Initializes the median filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public static Morphology Median(int radius)
        {
            int threshold = radius / 2;
            return new Morphology(radius, radius, threshold);
        }
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public static Morphology Erosion(int radius)
        {
            return new Morphology(radius, radius, 0);
        }
        /// <summary>
        /// Initializes the dilatation filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public static Morphology Dilatation(int radius)
        {
            return new Morphology(radius, radius, radius - 1);
        }
        #endregion
    }
}
