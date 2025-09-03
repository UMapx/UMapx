using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines Frei-Chen convolution filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Frei-Chen_operator
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FreiChen : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private FreiChenMode mode;
        private static readonly float[][] masks = CreateMasks();
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes Frei-Chen filter.
        /// </summary>
        /// <param name="mode">Filter mode</param>
        public FreiChen(FreiChenMode mode = FreiChenMode.Edge)
        {
            this.mode = mode;
        }
        /// <summary>
        /// Gets or sets filter mode.
        /// </summary>
        public FreiChenMode Mode
        {
            get => mode;
            set => mode = value;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            int width = bmData.Width;
            int height = bmData.Height;
            int stride = bmData.Stride;

            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();

            // Preserve alpha, copy layout once
            Buffer.MemoryCopy(src, dst, (long)height * stride, (long)height * stride);

            Parallel.For(0, height, y =>
            {
                float[] region = new float[9];
                float[] g = new float[9];
                int ys = y * stride;

                for (int x = 0; x < width; x++)
                {
                    // 3x3 neighborhood with replicate padding
                    int k = 0;
                    for (int di = -1; di <= 1; di++)
                    {
                        int yi = y + di;
                        if (yi < 0) yi = 0; else if (yi >= height) yi = height - 1;
                        int yoff = yi * stride;

                        for (int dj = -1; dj <= 1; dj++)
                        {
                            int xj = x + dj;
                            if (xj < 0) xj = 0; else if (xj >= width) xj = width - 1;

                            byte* p = &src[yoff + (xj * 4)];
                            region[k++] = 0.299f * p[2] + 0.587f * p[1] + 0.114f * p[0];
                        }
                    }

                    // project onto 9 Frei–Chen masks
                    for (int m = 0; m < 9; m++)
                    {
                        float s = 0f;
                        float[] mask = masks[m];
                        for (int t = 0; t < 9; t++) s += mask[t] * region[t];
                        g[m] = s;
                    }

                    float v;
                    switch (mode)
                    {
                        case FreiChenMode.Gradient:
                            v = MathF.Sqrt(g[0] * g[0] + g[1] * g[1]);
                            break;

                        case FreiChenMode.Edge:
                            v = MathF.Sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2] + g[3] * g[3]);
                            break;

                        case FreiChenMode.Line:
                            v = MathF.Sqrt(g[4] * g[4] + g[5] * g[5] + g[6] * g[6] + g[7] * g[7]);
                            break;

                        default: // FreiChenComponent.Average
                                 // g[8] is 3 * mean (because F9 is normalized with 1/3 entries)
                            v = g[8] * (1f / 3f);
                            break;
                    }

                    byte val = MathF.Byte(v);
                    int pos = ys + (x * 4);
                    dst[pos + 0] = val; // B
                    dst[pos + 1] = val; // G
                    dst[pos + 2] = val; // R
                                        // A preserved
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
            var Src = (Bitmap)Data.Clone();
            Apply(Data, Src);
            Src.Dispose();
        }
        #endregion

        #region Masks
        /// <summary>
        /// Creates masks.
        /// </summary>
        /// <returns>Jagged array</returns>
        private static float[][] CreateMasks()
        {
            float s2 = MathF.Sqrt2;
            float a = 0.5f;
            float b = 1.0f / 6.0f;
            return new float[][]
            {
                new float[] { a, a*s2, a, 0,0,0, -a,-a*s2,-a },
                new float[] { a,0,-a, a*s2,0,-a*s2, a,0,-a },
                new float[] { 0,-a,a*s2, a,0,-a, -a*s2,a,0 },
                new float[] { a*s2,-a,0, -a,0,a, 0,a,-a*s2 },
                new float[] { -b,-b*s2,-b, 2*b,2*b*s2,2*b, -b,-b*s2,-b },
                new float[] { -b,2*b,-b, -b*s2,2*b*s2,-b*s2, -b,2*b,-b },
                new float[] { 2*b,-b,-b*s2, -b,2*b,-b, -b*s2,-b,2*b },
                new float[] { -b*s2,-b,2*b, -b,2*b,-b, 2*b,-b,-b*s2 },
                new float[] { 1f/3f,1f/3f,1f/3f,1f/3f,1f/3f,1f/3f,1f/3f,1f/3f,1f/3f }
            };
        }
        #endregion
    }
}