using System;
using System.Drawing;
using System.Drawing.Imaging;
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
            Matrix = Core.Matrice.One(3, 3); Offset = 0; Bilateral = false;
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
        /// Prepares kernel dimensions and internal buffers.
        /// </summary>
        /// <param name="m">Convolution kernel</param>
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
            if (bmData.Width != bmSrc.Width || bmData.Height != bmSrc.Height)
                throw new ArgumentException("Bitmap sizes must match");

            if (bmData.PixelFormat != PixelFormat.Format32bppArgb || bmSrc.PixelFormat != PixelFormat.Format32bppArgb)
                throw new NotSupportedException("Only support Format32bppArgb pixelFormat");

            #region Data
            if (this.l0 != this.l1 && this.bilateral == true)
                throw new ArgumentException("Matrix must be squared");

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

        #region Sobel's gradient components
        /// <summary>
        /// Gets the value of the gradient operator.
        /// </summary>
        /// <param name="Gx">Gradient X</param>
        /// <param name="Gy">Gradient Y</param>
        /// <returns>Value</returns>
        public static float G(float Gx, float Gy)
        {
            return Maths.Sqrt(Gx * Gx + Gy * Gy);
        }
        /// <summary>
        /// Gets the angle of the gradient operator.
        /// </summary>
        /// <param name="Gx">Gradient X</param>
        /// <param name="Gy">Gradient Y</param>
        /// <returns>Value</returns>
        public static float Tetta(float Gx, float Gy)
        {
            return Maths.Atan(Gx / Gy);
        }
        #endregion

        #region Radius matrix
        /// <summary>
        /// Implements the construction of the Gaussian blur filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigmaX">Standard deviation X (>0)</param>
        /// <param name="sigmaY">Standard deviation Y (>0)</param>
        /// <returns>Matrix</returns>
        public static Convolution Gaussian(int m, int l, float sigmaY, float sigmaX)
        {
            return new Convolution(Operator.Gaussian(m, l, sigmaY, sigmaX));
        }
        /// <summary>
        /// Implements the construction of the "unsharp masking" filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigmaX">Standard deviation X (>0)</param>
        /// <param name="sigmaY">Standard deviation Y (>0)</param>
        /// <returns>Matrix</returns>
        public static Convolution Unsharp(int m, int l, float sigmaY, float sigmaX)
        {
            return new Convolution(Operator.Unsharp(m, l, sigmaY, sigmaX));
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
            return new Convolution(Operator.HighPass(m, l, boost));
        }
        /// <summary>
        /// Implements the construction of the low-pass filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Convolution LowPass(int m, int l)
        {
            return new Convolution(Operator.LowPass(m, l));
        }
        /// <summary>
        /// Implements the construction of the emboss filter.
        /// </summary>
        /// <param name="radius">Size</param>
        /// <returns>Matrix</returns>
        public static Convolution Emboss(int radius)
        {
            return new Convolution(Operator.Emboss(radius));
        }
        /// <summary>
        /// Implements the motion blur filter.
        /// </summary>
        /// <param name="radius">Size</param>
        /// <param name="angle">Angle in degrees</param>
        /// <param name="blur">Edge blur factor [0, 1]</param>
        public static Convolution MotionBlur(int radius, float angle, float blur)
        {
            return new Convolution(Operator.MotionBlur(radius, angle, blur));
        }
        #endregion

        #region Fixed radius matrix
        /// <summary>
        /// Implements the construction of the Roberts operator [2 x 2].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Roberts()
        {
            return new Convolution(Operator.Roberts());
        }
        /// <summary>
        /// Implements the construction of the Prewitt operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Prewitt()
        {
            return new Convolution(Operator.Prewitt());
        }
        /// <summary>
        /// Implements the construction of the Sobel operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Sobel()
        {
            return new Convolution(Operator.Sobel());
        }
        /// <summary>
        /// Implements the construction of the Scharr operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Scharr()
        {
            return new Convolution(Operator.Scharr());
        }
        /// <summary>
        /// Implements the construction of the Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Laplacian()
        {
            return new Convolution(Operator.Laplacian());
        }
        /// <summary>
        /// Implements the construction of the diagonal Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution LaplacianDiagonal()
        {
            return new Convolution(Operator.LaplacianDiagonal());
        }
        /// <summary>
        /// Implements the construction of the inverted Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution LaplacianInvert()
        {
            return new Convolution(Operator.LaplacianInvert());
        }
        #endregion

        #region Fixed radius compass matrix
        /// <summary>
        /// Implements the construction of the Kirsch operator [3 x 3].
        /// </summary>
        /// <param name="direction">Gradient direction</param>
        /// <returns>Matrix</returns>
        public static Convolution Kirsch(Gradient direction)
        {
            return new Convolution(Operator.Kirsch(direction));
        }
        /// <summary>
        /// Implements the construction of the Roberts operator [3 x 3]. [2 x 2].
        /// </summary>
        /// <param name="direction">Gradient direction</param>
        /// <returns>Matrix</returns>
        public static Convolution Roberts(Gradient direction)
        {
            return new Convolution(Operator.Roberts(direction));
        }
        #endregion
    }
}
