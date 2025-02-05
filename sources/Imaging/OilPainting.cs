using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the oil filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.codeproject.com/articles/471994/oilpainteffect
    /// </remarks>
    /// </summary>
    [Serializable]
    public class OilPainting : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private double depth;
        private int radius0;
        private int radius1;
        private int l0;
        private int l1;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the oil filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="depth">Value [0, 1]</param>
        public OilPainting(int radius = 3, double depth = 1.0)
        {
            Size = new SizeInt(radius, radius);
            Depth = depth;
        }
        /// <summary>
        /// Initializes the oil filter.
        /// </summary>
        /// <param name="height">Filter height</param>
        /// <param name="width">Filter width</param>
        /// <param name="depth">Value [0, 1]</param>
        public OilPainting(int width, int height, double depth = 1.0)
        {
            Size = new SizeInt(width, height);
            Depth = depth;
        }
        /// <summary>
        /// Initializes the oil filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="depth">Value [0, 1]</param>
        public OilPainting(SizeInt size, double depth = 1.0)
        {
            Size = size;
            Depth = depth;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(this.l0, this.l1);
            }
            set
            {
                this.l0 = value.Width;
                this.l1 = value.Height;
                this.radius0 = l0 >> 1;
                this.radius1 = l1 >> 1;
            }
        }
        /// <summary>
        /// Gets or sets the depth value [0, 1].
        /// </summary>
        public double Depth
        {
            get
            {
                return this.depth;
            }
            set
            {
                this.depth = value;
            }
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
            double strenght = this.depth * 255.0;
            double strenghtGlobal = this.depth / 3.0;
            #endregion

            Parallel.For(0, height, y =>
            {
                int[] Red = new int[256];
                int[] Green = new int[256];
                int[] Blue = new int[256];
                int[] Intensity = new int[256];

                int red, green, blue, intensity, max = 0, index = 0;
                int x, i, j, ir, jr, yr, xr, n;
                int ystride, v;
                byte* p;

                yr = y - radius0;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x - radius1;
                    v = ystride + x * 4;
                    red = green = blue = 0;

                    #region Convolution filtering
                    for (i = 0; i < l0; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        for (j = 0; j < l1; j++)
                        {
                            jr = xr + j;

                            if (jr < 0) continue; if (jr >= width) break;

                            p = &src[ir * stride + jr * 4];
                            red = p[2]; green = p[1]; blue = p[0];
                            intensity = (int)((red + green + blue) * strenghtGlobal);

                            Red[intensity] += red;
                            Green[intensity] += green;
                            Blue[intensity] += blue;
                            Intensity[intensity]++;
                        }
                    }
                    #endregion

                    #region Frequent intensity
                    for (n = 0; n < strenght; n++)
                    {
                        if (Intensity[n] > max)
                        {
                            max = Intensity[n];
                            index = n;
                        }
                    }
                    #endregion

                    #region Result pixel calculation
                    if (max != 0)
                    {
                        blue = Blue[index] / max;
                        green = Green[index] / max;
                        red = Red[index] / max;
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = (byte)red;
                    dst[v + 1] = (byte)green;
                    dst[v] = (byte)blue;
                    #endregion

                    #region Clear data
                    Array.Clear(Red, 0, 256);
                    Array.Clear(Green, 0, 256);
                    Array.Clear(Blue, 0, 256);
                    Array.Clear(Intensity, 0, 256);
                    max = 0; index = 0;
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
