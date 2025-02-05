using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the convolution filter.
    /// </summary>
    [Serializable]
    public class Convolution : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private float[][] kernel;
        private float offset;
        private int radius0;
        private int radius1;
        private int l0;
        private int l1;
        private bool bilateral;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the convolution filter.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="offset">Offset</param>
        /// <param name="bilateral">Bilateral processing or not</param>
        public Convolution(float[,] m, float offset = 0, bool bilateral = false)
        {
            Matrix = m; Offset = offset; Bilateral = bilateral;
        }
        /// <summary>
        /// Initializes the convolution filter.
        /// </summary>
        public Convolution()
        {
            Matrix = Matrice.One(3, 3); Offset = 0; Bilateral = false;
        }
        /// <summary>
        /// Gets or sets the convolution matrix.
        /// </summary>
        public float[,] Matrix
        {
            get
            {
                return Jagged.FromJagged(this.kernel);
            }
            set
            {
                Data(value);
            }
        }
        /// <summary>
        /// Gets or sets the offset value.
        /// </summary>
        public float Offset
        {
            get
            {
                return this.offset;
            }
            set
            {
                this.offset = value;
            }
        }
        /// <summary>
        /// Bilateral processing or not.
        /// </summary>
        public bool Bilateral
        {
            get
            {
                return this.bilateral;
            }
            set
            {
                this.bilateral = value;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Matrix</param>
        private void Data(float[,] m)
        {
            this.l0 = m.GetLength(0);
            this.l1 = m.GetLength(1);
            this.radius0 = l0 >> 1;
            this.radius1 = l1 >> 1;
            this.kernel = Jagged.ToJagged(m);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            if (this.l0 != this.l1 && this.bilateral == true)
                throw new Exception("Matrix must be squared");

            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            #endregion

            if (!this.bilateral)
            {
                Parallel.For(0, height, y =>
                {
                    float red, green, blue, div, k;
                    int x, i, j, ir, jr, yr, xr;
                    int ystride, v;
                    byte* p;

                    yr = y - radius0;
                    ystride = y * stride;

                    for (x = 0; x < width; x++)
                    {
                        xr = x - radius1;
                        v = ystride + x * 4;

                        red = green = blue = div = 0;

                        #region Convolution filtering
                        for (i = 0; i < l0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < l1; j++)
                            {
                                jr = xr + j;

                                if (jr < 0) continue; if (jr >= width) break;

                                k = kernel[i][j];

                                if (k != 0)
                                {
                                    p = &src[ir * stride + jr * 4];
                                    div += k;
                                    red += k * p[2]; green += k * p[1]; blue += k * p[0];
                                }
                            }
                        }
                        #endregion

                        #region Divider and offset
                        if (div != 0)
                        {
                            red /= div;
                            green /= div;
                            blue /= div;
                        }
                        if (offset != 0)
                        {
                            red += offset;
                            green += offset;
                            blue += offset;
                        }
                        #endregion

                        #region Recording pixel
                        dst[v + 2] = Maths.Byte(red);
                        dst[v + 1] = Maths.Byte(green);
                        dst[v] = Maths.Byte(blue);
                        #endregion
                    }
                });
            }
            else
            {
                Parallel.For(0, height, y =>
                {
                    float red1, green1, blue1, div1, k1;
                    float red2, green2, blue2, div2, k2;
                    int x, i, j, ir, jr, yr, xr;
                    int ystride, v;
                    byte* p;

                    yr = y - radius0;
                    ystride = y * stride;

                    for (x = 0; x < width; x++)
                    {
                        xr = x - radius1;
                        v = ystride + x * 4;

                        red1 = green1 = blue1 = div1 = 0;
                        red2 = green2 = blue2 = div2 = 0;

                        #region Convolution filtering
                        for (i = 0; i < l0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < l1; j++)
                            {
                                jr = xr + j;

                                if (jr < 0) continue; if (jr >= width) break;

                                k1 = kernel[i][j]; k2 = kernel[j][i];

                                p = &src[ir * stride + jr * 4];

                                div1 += k1; div2 += k2;
                                red1 += k1 * p[2]; green1 += k1 * p[1]; blue1 += k1 * p[0];
                                red2 += k2 * p[2]; green2 += k2 * p[1]; blue2 += k2 * p[0];
                            }
                        }
                        #endregion

                        #region Divider and offset
                        if (div1 != 0)
                        {
                            red1 /= div1;
                            green1 /= div1;
                            blue1 /= div1;
                        }
                        if (div2 != 0)
                        {
                            red2 /= div2;
                            green2 /= div2;
                            blue2 /= div2;
                        }
                        if (offset != 0)
                        {
                            red1 += offset;
                            green1 += offset;
                            blue1 += offset;

                            red2 += offset;
                            green2 += offset;
                            blue2 += offset;
                        }
                        #endregion

                        #region Recording pixel
                        dst[v + 2] = Maths.Byte(G(red1, red2));
                        dst[v + 1] = Maths.Byte(G(green1, green2));
                        dst[v] = Maths.Byte(G(blue1, blue2));
                        #endregion
                    }
                });
            }

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

        #region Sobel's gradient components
        /// <summary>
        /// Gets the value of the gradient operator.
        /// </summary>
        /// <param name="Gx">Gradient X</param>
        /// <param name="Gy">Gradient Y</param>
        /// <returns>Double precision floating point number</returns>
        private static float G(float Gx, float Gy)
        {
            return (float)Math.Sqrt(Gx * Gx + Gy * Gy);
        }
        /// <summary>
        /// Gets the angle of the gradient operator.
        /// </summary>
        /// <param name="Gx">Gradient X</param>
        /// <param name="Gy">Gradient Y</param>
        /// <returns>Double precision floating point number</returns>
        private static float Tetta(float Gx, float Gy)
        {
            return (float)Math.Atan(Gx / Gy);
        }
        #endregion

        #region Radius matrix
        /// <summary>
        /// Implements the construction of the inverted Gausssian filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigma">Standard deviation (!=0)</param>
        /// <returns>Matrix</returns>
        public static Convolution LoGaussian(int m, int l, float sigma)
        {
            int r1 = m / 2;
            int r2 = l / 2;
            float[,] H = new float[m, l];
            float sigma2 = sigma * sigma;
            float f0 = -1.0f / (Maths.Pi * sigma2 * sigma2);
            float f1 = 2.0f * sigma2;
            float kernel;
            int i, j, x, y;

            for (y = -r1, i = 0; i < m; y++, i++)
            {
                for (x = -r2, j = 0; j < l; x++, j++)
                {
                    kernel = (x * x + y * y) / f1;
                    H[i, j] = f0 * (1.0f - kernel) * Maths.Exp(-kernel);
                }
            }
            return new Convolution(H);
        }
        /// <summary>
        /// Implements the construction of the Gaussian blur filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigma">Standard deviation (!=0)</param>
        /// <returns>Matrix</returns>
        public static Convolution Gaussian(int m, int l, float sigma)
        {
            int r1 = m / 2;
            int r2 = l / 2;
            float[,] H = new float[m, l];
            int i, j, x, y;

            for (y = -r1, i = 0; i < m; y++, i++)
            {
                for (x = -r2, j = 0; j < l; x++, j++)
                {
                    H[i, j] = Kernel.Gaussian(x, sigma) * Kernel.Gaussian(y, sigma);
                }
            }
            return new Convolution(H);
        }
        /// <summary>
        /// Implements the construction of the "unsharp masking" filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigma">Standard deviation (!=0)</param>
        /// <returns>Matrix</returns>
        public static Convolution Unsharp(int m, int l, float sigma)
        {
            float[,] G = new float[m, l];
            int i, j, x, y;
            int r1 = m / 2;
            int r2 = l / 2;

            for (y = -r1, i = 0; i < m; y++, i++)
            {
                for (x = -r2, j = 0; j < l; x++, j++)
                {
                    G[i, j] = Kernel.Gaussian(x, sigma) * Kernel.Gaussian(y, sigma);
                }
            }

            float[,] invG = new float[m, l];
            float max = G[0, 0];
            float summary = 0;
            float v, iv;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    v = G[i, j] / max;
                    iv = v > byte.MaxValue ? byte.MaxValue : (int)v;
                    invG[i, j] = -iv;
                    summary += iv;
                }
            }

            invG[m / 2, l / 2] = 2 * summary - invG[m / 2, l / 2];
            return new Convolution(invG);
        }
        /// <summary>
        /// Implements the construction of the high-pass filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="boost">Boost</param>
        /// <returns>Matrix</returns>
        public static Convolution HighPass(int m, int l, float boost)
        {
            int r1 = m / 2;
            int r2 = l / 2;
            float[,] H = new float[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = -1;
                }
            }
            H[r1, r2] = boost;
            return new Convolution(H);
        }
        /// <summary>
        /// Implements the construction of the low-pass filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Convolution LowPass(int m, int l)
        {
            return Convolution.HighPass(m, l, 1);
        }
        /// <summary>
        /// Implements the construction of the emboss filter.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static Convolution Emboss(int n)
        {
            float[,] H = new float[n, n];
            int r = n - 1, r2 = r / 2;

            H[0, 0] = -2; H[0, r2] = -1; //         0;
            H[r2, 0] = -1; H[r2, r2] = 1; H[r2, r] = 1;
            H[r, r2] = 1; H[r, r] = 2; //         0;

            return new Convolution(H);
        }
        #endregion

        #region Fixed radius matrix
        /// <summary>
        /// Implements the construction of the Roberts operator [2 x 2].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Roberts()
        {
            return new Convolution(new float[2, 2] { { 1, 0 }, { 0, -1 } });
        }
        /// <summary>
        /// Implements the construction of the Prewitt operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Prewitt()
        {
            return new Convolution(new float[3, 3] { { -1, -1, -1 }, { 0, 0, 0 }, { 1, 1, 1 } });
        }
        /// <summary>
        /// Implements the construction of the Sobel operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Sobel()
        {
            return new Convolution(new float[3, 3] { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } });
        }
        /// <summary>
        /// Implements the construction of the Scharr operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Scharr()
        {
            return new Convolution(new float[3, 3] { { 3, 10, 3 }, { 0, 0, 0 }, { -3, -10, -3 } });
        }
        /// <summary>
        /// Implements the construction of the Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Laplacian()
        {
            return new Convolution(new float[3, 3] { { 0, 1, 0 }, { 1, -4, 1 }, { 0, 1, 0 } });
        }
        /// <summary>
        /// Implements the construction of the diagonal Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution LaplacianDiagonal()
        {
            return new Convolution(new float[3, 3] { { 1, 1, 1 }, { 1, -8, 1 }, { 1, 1, 1 } });
        }
        /// <summary>
        /// Implements the construction of the inverted Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution LaplacianInvert()
        {
            return new Convolution(new float[3, 3] { { -1, 0, -1 }, { 0, 4, 0 }, { -1, 0, -1 } });
        }
        #endregion

        #region Fixed radius compass matrix
        /// <summary>
        /// Implements the construction of the Kirsh operator [3 x 3].
        /// </summary>
        /// <param name="direction">Gradient direction</param>
        /// <returns>Matrix</returns>
        public static Convolution Kirsh(Gradient direction)
        {
            float[,] H = new float[3, 3];

            if (direction == Gradient.North)
            {
                H[0, 0] = 5; H[0, 1] = 5; H[0, 2] = 5;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = -3;
            }
            else if (direction == Gradient.NorthWest)
            {
                H[0, 0] = 5; H[0, 1] = 5; H[0, 2] = -3;
                H[1, 0] = 5; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = -3;
            }
            else if (direction == Gradient.West)
            {
                H[0, 0] = 5; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = 5; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = 5; H[2, 1] = -3; H[2, 2] = -3;
            }
            else if (direction == Gradient.SouthWest)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = 5; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = 5; H[2, 1] = 5; H[2, 2] = -3;
            }
            else if (direction == Gradient.South)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = 5; H[2, 1] = 5; H[2, 2] = 5;
            }
            else if (direction == Gradient.SouthEast)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = 5;
                H[2, 0] = -3; H[2, 1] = 5; H[2, 2] = 5;
            }
            else if (direction == Gradient.East)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = 5;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = 5;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = 5;
            }
            else if (direction == Gradient.NorthEast)
            {
                H[0, 0] = -3; H[0, 1] = 5; H[0, 2] = 5;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = 5;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = -3;
            }

            return new Convolution(H);
        }
        /// <summary>
        /// Implements the construction of the Roberts operator [3 x 3]. [2 x 2].
        /// </summary>
        /// <param name="direction">Gradient direction</param>
        /// <returns>Matrix</returns>
        public static Convolution Roberts(Gradient direction)
        {
            float[,] H = new float[2, 2];

            if (direction == Gradient.North)
            {
                H[0, 0] = 0; H[0, 1] = 1;
                H[1, 0] = -1; H[1, 1] = 0;
            }
            else if (direction == Gradient.West)
            {
                H[0, 0] = 1; H[0, 1] = 0;
                H[1, 0] = 0; H[1, 1] = -1;
            }
            else if (direction == Gradient.South)
            {
                H[0, 0] = 0; H[0, 1] = -1;
                H[1, 0] = 1; H[1, 1] = 0;
            }
            else // Properties.Gradient.East
            {
                H[0, 0] = -1; H[0, 1] = 0;
                H[1, 0] = 0; H[1, 1] = 1;
            }

            return new Convolution(H);
        }
        #endregion
    }
}
