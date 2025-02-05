using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Uses to work with bitmap matrices.
    /// </summary>
    public static class BitmapMatrix
    {
        #region RGB
        /// <summary>
        /// Converts a Bitmap to an RGB structure with or without alpha-channel.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>RGBA structure array</returns>
        public static float[][,] ToRGB(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            float[][,] rgb = BitmapMatrix.ToRGB(bmData, alpha);
            BitmapFormat.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an RGB structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>RGBA structure array</returns>
        public unsafe static float[][,] ToRGB(BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            float[,] x = new float[height, width];
            float[,] y = new float[height, width];
            float[,] z = new float[height, width];

            // with alpha channel
            if (alpha)
            {
                float[,] a = new float[height, width];

                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // color
                        x[j, i] = p[k + 0] / 255.0f;
                        y[j, i] = p[k + 1] / 255.0f;
                        z[j, i] = p[k + 2] / 255.0f;
                        a[j, i] = p[k + 3] / 255.0f;
                    }
                });

                return new float[][,] { x, y, z, a };
            }

            // without alpha channel
            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * 4;

                    // color
                    x[j, i] = p[k + 0] / 255.0f;
                    y[j, i] = p[k + 1] / 255.0f;
                    z[j, i] = p[k + 2] / 255.0f;
                }
            });

            return new float[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an RGB structure to a color image.
        /// </summary>
        /// <param name="array">RGBA structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromRGB(this float[][,] array)
        {
            // matrices
            float[,] x = array[0];
            float[,] y = array[1];
            float[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapFormat.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                float[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0f);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0f);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0f);
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0f);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0f);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0f);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0f);
                        p[k + 3] = (byte)255;
                    }
                });
            }


            BitmapFormat.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an RGB structure to a color image.
        /// </summary>
        /// <param name="array">RGBA structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromRGB(this float[][,] array, BitmapData bmData)
        {
            // matrices
            float[,] x = array[0];
            float[,] y = array[1];
            float[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                float[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0f);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0f);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0f);
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0f);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0f);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0f);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0f);
                        p[k + 3] = (byte)255;
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Converts an RGB structure to a color image.
        /// </summary>
        /// <param name="array">RGBA structure array</param>
        /// <param name="Data">Bitmap</param>
        public static void FromRGB(this float[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            FromRGB(array, bmData);
            BitmapFormat.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region HSB
        /// <summary>
        /// Converts a Bitmap to an HSB structure with or without alpha-channel.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSB structure array</returns>
        public static float[][,] ToHSB(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            float[][,] rgb = BitmapMatrix.ToHSB(bmData, alpha);
            BitmapFormat.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an HSB structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSB structure array</returns>
        public unsafe static float[][,] ToHSB(this BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            float[,] x = new float[height, width];
            float[,] y = new float[height, width];
            float[,] z = new float[height, width];

            // with alpha channel
            if (alpha)
            {
                float[,] a = new float[height, width];

                Parallel.For(0, height, j =>
                {
                    HSB mdl;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        mdl = HSB.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                        x[j, i] = mdl.Hue;
                        y[j, i] = mdl.Saturation;
                        z[j, i] = mdl.Brightness;
                        a[j, i] = p[k + 3] / 255.0f;
                    }
                });

                return new float[][,] { x, y, z, a };
            }

            // without alpha channel
            Parallel.For(0, height, j =>
            {
                HSB mdl;

                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * 4;

                    mdl = HSB.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    x[j, i] = mdl.Hue;
                    y[j, i] = mdl.Saturation;
                    z[j, i] = mdl.Brightness;
                }
            });

            return new float[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an HSB structure to a color image.
        /// </summary>
        /// <param name="array">HSB structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromHSB(this float[][,] array)
        {
            // matrices
            float[,] x = array[0];
            float[,] y = array[1];
            float[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapFormat.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                float[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSB(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0f);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSB(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            BitmapFormat.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an HSB structure to a color image.
        /// </summary>
        /// <param name="array">HSB structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromHSB(this float[][,] array, BitmapData bmData)
        {
            // matrices
            float[,] x = array[0];
            float[,] y = array[1];
            float[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                float[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSB(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0f);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSB(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Converts an HSB structure to a color image.
        /// </summary>
        /// <param name="array">HSB structure array</param>
        /// <param name="Data">Bitmap</param>
        public static void FromHSB(this float[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            FromHSB(array, bmData);
            BitmapFormat.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region HSL
        /// <summary>
        /// Converts a Bitmap to an HSL structure with or without alpha-channel.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSL structure array</returns>
        public static float[][,] ToHSL(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            float[][,] rgb = BitmapMatrix.ToHSL(bmData, alpha);
            BitmapFormat.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an HSL structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSL structure array</returns>
        public unsafe static float[][,] ToHSL(this BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            float[,] x = new float[height, width];
            float[,] y = new float[height, width];
            float[,] z = new float[height, width];

            // with alpha channel
            if (alpha)
            {
                float[,] a = new float[height, width];

                Parallel.For(0, height, j =>
                {
                    HSL mdl;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        mdl = HSL.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                        x[j, i] = mdl.Hue;
                        y[j, i] = mdl.Saturation;
                        z[j, i] = mdl.Lightness;
                        a[j, i] = p[k + 3] / 255.0f;
                    }
                });

                return new float[][,] { x, y, z, a };
            }

            // without alpha channel
            Parallel.For(0, height, j =>
            {
                HSL mdl;

                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * 4;

                    mdl = HSL.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    x[j, i] = mdl.Hue;
                    y[j, i] = mdl.Saturation;
                    z[j, i] = mdl.Lightness;
                }
            });

            return new float[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an HSL structure to a color image.
        /// </summary>
        /// <param name="array">HSL structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromHSL(this float[][,] array)
        {
            // matrices
            float[,] x = array[0];
            float[,] y = array[1];
            float[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapFormat.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                float[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSL(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0f);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSL(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            BitmapFormat.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an HSL structure to a color image.
        /// </summary>
        /// <param name="array">HSL structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromHSL(this float[][,] array, BitmapData bmData)
        {
            // matrices
            float[,] x = array[0];
            float[,] y = array[1];
            float[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                float[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSL(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0f);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSL(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Converts an HSL structure to a color image.
        /// </summary>
        /// <param name="array">HSL structure array</param>
        /// <param name="Data">Bitmap</param>
        public static void FromHSL(this float[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            FromHSL(array, bmData);
            BitmapFormat.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region YCbCr
        /// <summary>
        /// Converts a Bitmap to an YCbCr structure with or without alpha-channel.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>YCbCr structure array</returns>
        public static float[][,] ToYCbCr(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            float[][,] rgb = BitmapMatrix.ToYCbCr(bmData, alpha);
            BitmapFormat.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an YCbCr structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>YCbCr structure array</returns>
        public unsafe static float[][,] ToYCbCr(this BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            float[,] x = new float[height, width];
            float[,] y = new float[height, width];
            float[,] z = new float[height, width];

            // with alpha channel
            if (alpha)
            {
                float[,] a = new float[height, width];

                Parallel.For(0, height, j =>
                {
                    YCbCr mdl;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        mdl = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                        x[j, i] = mdl.Y;
                        y[j, i] = mdl.Cb;
                        z[j, i] = mdl.Cr;
                        a[j, i] = p[k + 3] / 255.0f;
                    }
                });

                return new float[][,] { x, y, z, a };
            }

            // without alpha channel
            Parallel.For(0, height, j =>
            {
                YCbCr mdl;

                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * 4;

                    mdl = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    x[j, i] = mdl.Y;
                    y[j, i] = mdl.Cb;
                    z[j, i] = mdl.Cr;
                }
            });

            return new float[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an YCbCr structure to a color image.
        /// </summary>
        /// <param name="array">YCbCr structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromYCbCr(this float[][,] array)
        {
            // matrices
            float[,] x = array[0];
            float[,] y = array[1];
            float[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapFormat.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                float[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new YCbCr(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0f);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new YCbCr(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            BitmapFormat.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an YCbCr structure to a color image.
        /// </summary>
        /// <param name="array">YCbCr structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromYCbCr(this float[][,] array, BitmapData bmData)
        {
            // matrices
            float[,] x = array[0];
            float[,] y = array[1];
            float[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                float[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new YCbCr(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0f);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new YCbCr(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Converts an YCbCr structure to a color image.
        /// </summary>
        /// <param name="array">YCbCr structure array</param>
        /// <param name="Data">Bitmap</param>
        public static void FromYCbCr(this float[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            FromYCbCr(array, bmData);
            BitmapFormat.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region Grayscale
        /// <summary>
        /// Converts Bitmap to averaged channel value matrix.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Matrix</returns>
        public static float[,] ToGrayscale(this Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            float[,] rgb = ToGrayscale(bmData);
            BitmapFormat.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts Bitmap to averaged channel value matrix.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Matrix</returns>
        public unsafe static float[,] ToGrayscale(this BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float[,] rgb = new float[height, width];
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    rgb[j, i] = RGB.Average(p[k + 2], p[k + 1], p[k]) / 255.0f;
                }
            });

            return rgb;
        }
        /// <summary>
        /// Converts a matrix of channel values to a monochrome Bitmap.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromGrayscale(this float[,] m)
        {
            int width = m.GetLength(1), height = m.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapFormat.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(m[j, i] * 255.0f);
                    p[k + 3] = 255;
                }
            });

            BitmapFormat.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts a matrix of channel values to a monochrome Bitmap.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromGrayscale(this float[,] m, BitmapData bmData)
        {
            int width = m.GetLength(1), height = m.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(m[j, i] * 255.0f);
                    p[k + 3] = 255;
                }
            });

            return;
        }
        /// <summary>
        /// Converts a matrix of channel values to a monochrome Bitmap.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="Data">Bitmap</param>
        public static void FromGrayscale(this float[,] m, Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            FromGrayscale(m, bmData);
            BitmapFormat.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
}
