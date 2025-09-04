using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Kuwahara filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Kuwahara_filter
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Kuwahara : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int rw;
        private int rh;
        private int windows = 4;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the Kuwahara filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="windows">Number of windows (4 or 8)</param>
        public Kuwahara(int radius = 4, int windows = 4)
        {
            Size = new SizeInt(radius, radius);
            Windows = windows;
        }
        /// <summary>
        /// Initializes the Kuwahara filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="windows">Number of windows (4 or 8)</param>
        public Kuwahara(int width, int height, int windows = 4)
        {
            Size = new SizeInt(width, height);
            Windows = windows;
        }
        /// <summary>
        /// Initializes the Kuwahara filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="windows">Number of windows (4 or 8)</param>
        public Kuwahara(SizeInt size, int windows = 4)
        {
            Size = size;
            Windows = windows;
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
        /// Gets or sets the number of windows (4 or 8).
        /// </summary>
        public int Windows
        {
            get { return windows; }
            set { windows = Math.Max(4, Math.Min(8, value)); }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int rx = rw >> 1;
            int ry = rh >> 1;
            int winCount = windows < 8 ? 4 : 8;

            int[,] rects = winCount == 4 ? new int[,] {
                {-rx, -ry, 0, 0},
                {0, -ry, rx, 0},
                {-rx, 0, 0, ry},
                {0, 0, rx, ry}
            } : new int[,] {
                {-rx, -ry, 0, 0},
                {0, -ry, rx, 0},
                {-rx, 0, 0, ry},
                {0, 0, rx, ry},
                {-rx, -ry, rx, 0},
                {-rx, 0, rx, ry},
                {-rx, -ry, 0, ry},
                {0, -ry, rx, ry}
            };

            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    double minVar = double.MaxValue;
                    double outR = 0, outG = 0, outB = 0;

                    for (int k = 0; k < winCount; k++)
                    {
                        int x1 = Math.Max(0, x + rects[k, 0]);
                        int y1 = Math.Max(0, y + rects[k, 1]);
                        int x2 = Math.Min(width - 1, x + rects[k, 2]);
                        int y2 = Math.Min(height - 1, y + rects[k, 3]);

                        double sumR = 0, sumG = 0, sumB = 0, sumI = 0, sumI2 = 0;
                        int count = 0;

                        for (int yy = y1; yy <= y2; yy++)
                        {
                            byte* p = src + yy * stride + x1 * 4;
                            for (int xx = x1; xx <= x2; xx++, p += 4)
                            {
                                double r = p[2];
                                double g = p[1];
                                double b = p[0];
                                double i = (r + g + b) / 3.0;

                                sumR += r;
                                sumG += g;
                                sumB += b;
                                sumI += i;
                                sumI2 += i * i;
                                count++;
                            }
                        }

                        double meanI = sumI / count;
                        double variance = sumI2 / count - meanI * meanI;

                        if (variance < minVar)
                        {
                            minVar = variance;
                            outR = sumR / count;
                            outG = sumG / count;
                            outB = sumB / count;
                        }
                    }

                    int q = y * stride + x * 4;
                    dst[q + 2] = Maths.Byte((float)outR);
                    dst[q + 1] = Maths.Byte((float)outG);
                    dst[q] = Maths.Byte((float)outB);
                }
            });
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
        }
        #endregion
    }
}