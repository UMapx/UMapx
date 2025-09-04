using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Wiener filter for local noise reduction.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wiener_filter
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Wiener : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int rw;
        private int rh;
        private float noise;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the Wiener filter with square window.
        /// </summary>
        /// <param name="radius">Window radius</param>
        /// <param name="noise">Noise variance</param>
        public Wiener(int radius = 3, float noise = 400f)
        {
            Size = new SizeInt(radius, radius);
            Noise = noise;
        }
        /// <summary>
        /// Initializes the Wiener filter.
        /// </summary>
        /// <param name="width">Window width</param>
        /// <param name="height">Window height</param>
        /// <param name="noise">Noise variance</param>
        public Wiener(int width, int height, float noise = 400f)
        {
            Size = new SizeInt(width, height);
            Noise = noise;
        }
        /// <summary>
        /// Initializes the Wiener filter.
        /// </summary>
        /// <param name="size">Window size</param>
        /// <param name="noise">Noise variance</param>
        public Wiener(SizeInt size, float noise = 400f)
        {
            Size = size;
            Noise = noise;
        }
        /// <summary>
        /// Gets or sets the window size.
        /// </summary>
        public SizeInt Size
        {
            get { return new SizeInt(rw, rh); }
            set { rw = value.Width; rh = value.Height; }
        }
        /// <summary>
        /// Gets or sets the noise variance.
        /// </summary>
        public float Noise
        {
            get { return noise; }
            set { noise = value; }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            ApplyInternal(bmData, bmSrc);
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
            Bitmap src = (Bitmap)current.Clone();
            BitmapData bmSrc = BitmapFormat.Lock32bpp(src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(src, bmSrc);
            src.Dispose();
            current.Dispose();
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            var src = (Bitmap)Data.Clone();
            Apply(Data, src);
            src.Dispose();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Internal apply method.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyInternal(BitmapData bmData, BitmapData bmSrc)
        {
            int width = bmData.Width;
            int height = bmData.Height;
            int stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();

            int halfW = rw / 2;
            int halfH = rh / 2;
            int winW = 2 * halfW + 1;
            int winH = 2 * halfH + 1;
            double area = (double)winW * winH;

            float nv = Math.Max(0f, noise); // guard

            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    double sumR = 0, sumG = 0, sumB = 0;
                    double sumR2 = 0, sumG2 = 0, sumB2 = 0;

                    for (int j = -halfH; j <= halfH; j++)
                    {
                        int yy = MathsF.Range(y + j, 0, height - 1);
                        byte* pRow = src + yy * stride;

                        for (int i = -halfW; i <= halfW; i++)
                        {
                            int xx = MathsF.Range(x + i, 0, width - 1);
                            byte* p = pRow + xx * 4;

                            double b = p[0];
                            double g = p[1];
                            double r = p[2];

                            sumR += r; sumG += g; sumB += b;
                            sumR2 += r * r; sumG2 += g * g; sumB2 += b * b;
                        }
                    }

                    double meanR = sumR / area;
                    double meanG = sumG / area;
                    double meanB = sumB / area;

                    // numerically safe variance (can be ~-1e-12 from rounding)
                    double varR = Math.Max(0.0, sumR2 / area - meanR * meanR);
                    double varG = Math.Max(0.0, sumG2 / area - meanG * meanG);
                    double varB = Math.Max(0.0, sumB2 / area - meanB * meanB);

                    byte* pDst = dst + y * stride + x * 4;
                    byte* pSrc = src + y * stride + x * 4;

                    double r0 = pSrc[2];
                    double g0 = pSrc[1];
                    double b0 = pSrc[0];

                    // Wiener shrinkage: gain in [0,1]
                    double gainR = (varR > 0.0) ? Math.Max(0.0, varR - nv) / varR : 0.0;
                    double gainG = (varG > 0.0) ? Math.Max(0.0, varG - nv) / varG : 0.0;
                    double gainB = (varB > 0.0) ? Math.Max(0.0, varB - nv) / varB : 0.0;

                    pDst[2] = MathsF.Byte((float)(meanR + gainR * (r0 - meanR)));
                    pDst[1] = MathsF.Byte((float)(meanG + gainG * (g0 - meanG)));
                    pDst[0] = MathsF.Byte((float)(meanB + gainB * (b0 - meanB)));
                    pDst[3] = pSrc[3];
                }
            });
        }
        #endregion
    }
}