using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines a normal bump filter (height → tangent-space normal map).
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Normal_mapping
    /// </remarks>
    [Serializable]
    public class NormalBump : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private float strength = 2.0f;
        private bool invertY = true;
        private bool useSobel = true;
        private bool wrapEdges = true;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the normal bump filter.
        /// </summary>
        /// <param name="strength">Strength (recommended 0.5..8.0)</param>
        public NormalBump(int strength = 0)
        {
            Strength = 1.0f + strength;
        }
        /// <summary>
        /// Gets or sets bump strength (recommended 0.5..8.0).
        /// </summary>
        public float Strength
        {
            get => strength;
            set => strength = Maths.Max(0f, value);
        }
        /// <summary>
        /// Invert the Y/green channel (true for DirectX, false for OpenGL).
        /// </summary>
        public bool InvertY
        {
            get => invertY;
            set => invertY = value;
        }
        /// <summary>
        /// Use Sobel 3x3 (true) or central differences (false).
        /// </summary>
        public bool UseSobel
        {
            get => useSobel;
            set => useSobel = value;
        }
        /// <summary>
        /// Wrap edges for tiling (true) or clamp (false).
        /// </summary>
        public bool WrapEdges
        {
            get => wrapEdges;
            set => wrapEdges = value;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            if (bmData.Width != bmSrc.Width || bmData.Height != bmSrc.Height)
                throw new ArgumentException("Bitmap sizes must match");

            if (bmData.PixelFormat != PixelFormat.Format32bppArgb || bmSrc.PixelFormat != PixelFormat.Format32bppArgb)
                throw new NotSupportedException("Only support Format32bppArgb pixelFormat");

            int w = bmData.Width, h = bmData.Height;
            int srcBpp = 4, dstBpp = 4;

            unsafe
            {
                byte* srcBase = (byte*)bmSrc.Scan0;
                byte* dstBase = (byte*)bmData.Scan0;

                // Sobel kernels
                int[,] Kx = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
                int[,] Ky = { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };

                (int x, int y) edge(int x, int y)
                {
                    if (wrapEdges)
                    {
                        if (x < 0) x += w; else if (x >= w) x -= w;
                        if (y < 0) y += h; else if (y >= h) y -= h;
                    }
                    else
                    {
                        if (x < 0) x = 0; else if (x >= w) x = w - 1;
                        if (y < 0) y = 0; else if (y >= h) y = h - 1;
                    }
                    return (x, y);
                }

                float sampleGray(int x, int y)
                {
                    var (ix, iy) = edge(x, y);
                    byte* p = srcBase + iy * bmSrc.Stride + ix * srcBpp;
                    float r = p[2] / 255f, g = p[1] / 255f, b = p[0] / 255f;
                    return 0.2126f * r + 0.7152f * g + 0.0722f * b;
                }

                Parallel.For(0, h, y =>
                {
                    byte* outRow = dstBase + y * bmData.Stride;

                    for (int x = 0; x < w; x++)
                    {
                        float dx, dy;

                        if (useSobel)
                        {
                            float gx = 0, gy = 0;
                            for (int j = -1; j <= 1; j++)
                                for (int i = -1; i <= 1; i++)
                                {
                                    float v = sampleGray(x + i, y + j);
                                    gx += v * Kx[j + 1, i + 1];
                                    gy += v * Ky[j + 1, i + 1];
                                }
                            dx = gx; dy = gy;
                        }
                        else
                        {
                            float l = sampleGray(x - 1, y);
                            float r = sampleGray(x + 1, y);
                            float u = sampleGray(x, y - 1);
                            float d = sampleGray(x, y + 1);
                            dx = (r - l) * 0.5f;
                            dy = (d - u) * 0.5f;
                        }

                        float nx = -dx * strength;
                        float ny = dy * strength;
                        if (invertY) ny = -ny;

                        float nz = 1f;
                        float inv = 1.0f / Maths.Sqrt(nx * nx + ny * ny + nz * nz);
                        nx *= inv; ny *= inv; nz *= inv;

                        byte R = ToByte(nx);
                        byte G = ToByte(ny);
                        byte B = ToByte(nz);

                        int o = x * dstBpp;
                        outRow[o + 0] = B;
                        outRow[o + 1] = G;
                        outRow[o + 2] = R;
                        outRow[o + 3] = 255;
                    }
                });
            }

            static byte ToByte(float v)
            {
                int t = (int)Math.Round((v * 0.5f + 0.5f) * 255f);
                if (t < 0) t = 0; else if (t > 255) t = 255;
                return (byte)t;
            }
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
            Bitmap Src = BitmapFormat.ToBitmap(bmData);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Src, bmSrc);
            Src.Dispose();
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
