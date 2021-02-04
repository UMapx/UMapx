using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Uses to work with bitmaps.
    /// </summary>
    public static class BitmapConverter
    {
        #region Bitmap convert components
        /// <summary>
        /// Converts Bitmap to icon file.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="size">Size</param>
        /// <returns>Icon</returns>
        public static Icon ToIco(this Bitmap b, int size)
        {
            using Bitmap bmp = new Bitmap(b, new Size(size, size));
            using MemoryStream pngstream = new MemoryStream();

            byte[] pngicon = new byte[] { 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            byte[] png;

            bmp.Save(pngstream, ImageFormat.Png);
            pngstream.Position = 0;
            png = pngstream.ToArray();

            if (size >= 256) size = 0;
            pngicon[6] = (byte)size;
            pngicon[7] = (byte)size;
            pngicon[14] = (byte)(png.Length & 255);
            pngicon[15] = (byte)(png.Length / 256);
            pngicon[18] = (byte)(pngicon.Length);

            MemoryStream icostream = new MemoryStream();
            icostream.Write(pngicon, 0, pngicon.Length);
            icostream.Write(png, 0, png.Length);
            icostream.Position = 0;

            return new Icon(icostream);
        }
        /// <summary>
        /// Converts Bitmap to JPEG format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToJpeg(this Bitmap b)
        {
            MemoryStream stream = new MemoryStream();
            b.Save(stream, ImageFormat.Jpeg);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Converts Bitmap to BMP format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToBmp(this Bitmap b)
        {
            MemoryStream stream = new MemoryStream();
            b.Save(stream, ImageFormat.Bmp);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Converts Bitmap to GIF format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToGif(this Bitmap b)
        {
            MemoryStream stream = new MemoryStream();
            b.Save(stream, ImageFormat.Gif);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Converts Bitmap to PNG format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToPng(this Bitmap b)
        {
            MemoryStream stream = new MemoryStream();
            b.Save(stream, ImageFormat.Png);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Converts Bitmap to TIFF format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToTiff(this Bitmap b)
        {
            MemoryStream stream = new MemoryStream();
            b.Save(stream, ImageFormat.Tiff);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Gets the Bitmap from the BitmapData.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Bitmap(this BitmapData bmData)
        {
            return new Bitmap(bmData.Width, bmData.Height, bmData.Stride, bmData.PixelFormat, bmData.Scan0);
        }
        /// <summary>
        /// Converts Bitmap to a specific format
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="pixelformat">Pixel format</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Bitmap(this Bitmap b, PixelFormat pixelformat)
        {
            return b.Clone(new Rectangle(0, 0, b.Width, b.Height), pixelformat);
        }
        #endregion

        #region BitmapData voids
        /// <summary>
        /// Blocks Bitmap in system memory.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap data</returns>
        public static BitmapData Lock32bpp(this Bitmap b)
        {
            return b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format32bppArgb);
        }
        /// <summary>
        /// Blocks Bitmap in system memory.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap data</returns>
        internal static BitmapData Lock8bpp(this Bitmap b)
        {
            return b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format8bppIndexed);
        }
        /// <summary>
        /// Unblocks Bitmap in system memory.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="bmData">Bitmap data</param>
        public static void Unlock(this Bitmap b, BitmapData bmData)
        {
            b.UnlockBits(bmData);
            return;
        }
        #endregion

        #region RGB
        /// <summary>
        /// Converts a Bitmap to an RGB structure with or without alpha-channel.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>RGBA structure array</returns>
        public static double[][,] ToRGB(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToRGB(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an RGB structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>RGBA structure array</returns>
        public unsafe static double[][,] ToRGB(BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            double[,] x = new double[height, width];
            double[,] y = new double[height, width];
            double[,] z = new double[height, width];

            // with alpha channel
            if (alpha)
            {
                double[,] a = new double[height, width];

                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // color
                        x[j, i] = p[k + 0] / 255.0;
                        y[j, i] = p[k + 1] / 255.0;
                        z[j, i] = p[k + 2] / 255.0;
                        a[j, i] = p[k + 3] / 255.0;
                    }
                });

                return new double[][,] { x, y, z, a };
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
                    x[j, i] = p[k + 0] / 255.0;
                    y[j, i] = p[k + 1] / 255.0;
                    z[j, i] = p[k + 2] / 255.0;
                }
            });

            return new double[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an RGB structure to a color image.
        /// </summary>
        /// <param name="array">RGBA structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromRGB(this double[][,] array)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0);
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
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
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0);
                        p[k + 3] = (byte)255;
                    }
                });
            }


            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an RGB structure to a color image.
        /// </summary>
        /// <param name="array">RGBA structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromRGB(this double[][,] array, BitmapData bmData)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0);
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
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
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0);
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
        public static void FromRGB(this double[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromRGB(array, bmData);
            BitmapConverter.Unlock(Data, bmData);
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
        public static double[][,] ToHSB(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToHSB(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an HSB structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSB structure array</returns>
        public unsafe static double[][,] ToHSB(this BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            double[,] x = new double[height, width];
            double[,] y = new double[height, width];
            double[,] z = new double[height, width];

            // with alpha channel
            if (alpha)
            {
                double[,] a = new double[height, width];

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
                        a[j, i] = p[k + 3] / 255.0;
                    }
                });

                return new double[][,] { x, y, z, a };
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

            return new double[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an HSB structure to a color image.
        /// </summary>
        /// <param name="array">HSB structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromHSB(this double[][,] array)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

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
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
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

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an HSB structure to a color image.
        /// </summary>
        /// <param name="array">HSB structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromHSB(this double[][,] array, BitmapData bmData)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

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
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
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
        public static void FromHSB(this double[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromHSB(array, bmData);
            BitmapConverter.Unlock(Data, bmData);
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
        public static double[][,] ToHSL(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToHSL(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an HSL structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSL structure array</returns>
        public unsafe static double[][,] ToHSL(this BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            double[,] x = new double[height, width];
            double[,] y = new double[height, width];
            double[,] z = new double[height, width];

            // with alpha channel
            if (alpha)
            {
                double[,] a = new double[height, width];

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
                        a[j, i] = p[k + 3] / 255.0;
                    }
                });

                return new double[][,] { x, y, z, a };
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

            return new double[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an HSL structure to a color image.
        /// </summary>
        /// <param name="array">HSL structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromHSL(this double[][,] array)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

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
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
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

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an HSL structure to a color image.
        /// </summary>
        /// <param name="array">HSL structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromHSL(this double[][,] array, BitmapData bmData)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

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
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
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
        public static void FromHSL(this double[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromHSL(array, bmData);
            BitmapConverter.Unlock(Data, bmData);
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
        public static double[][,] ToYCbCr(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToYCbCr(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an YCbCr structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>YCbCr structure array</returns>
        public unsafe static double[][,] ToYCbCr(this BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            double[,] x = new double[height, width];
            double[,] y = new double[height, width];
            double[,] z = new double[height, width];

            // with alpha channel
            if (alpha)
            {
                double[,] a = new double[height, width];

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
                        a[j, i] = p[k + 3] / 255.0;
                    }
                });

                return new double[][,] { x, y, z, a };
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

            return new double[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an YCbCr structure to a color image.
        /// </summary>
        /// <param name="array">YCbCr structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromYCbCr(this double[][,] array)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

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
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
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

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an YCbCr structure to a color image.
        /// </summary>
        /// <param name="array">YCbCr structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromYCbCr(this double[][,] array, BitmapData bmData)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

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
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
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
        public static void FromYCbCr(this double[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromYCbCr(array, bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region Grayscale
        /// <summary>
        /// Converts Bitmap to averaged channel value matrix.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Matrix</returns>
        public static double[,] ToGrayscale(this Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[,] rgb = ToGrayscale(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts Bitmap to averaged channel value matrix.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Matrix</returns>
        public unsafe static double[,] ToGrayscale(this BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double[,] rgb = new double[height, width];
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    rgb[j, i] = RGB.Average(p[k + 2], p[k + 1], p[k]) / 255.0;
                }
            });

            return rgb;
        }
        /// <summary>
        /// Converts a matrix of channel values to a monochrome Bitmap.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromGrayscale(this double[,] m)
        {
            int width = m.GetLength(1), height = m.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(m[j, i] * 255.0);
                    p[k + 3] = 255;
                }
            });

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts a matrix of channel values to a monochrome Bitmap.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromGrayscale(this double[,] m, BitmapData bmData)
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
                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(m[j, i] * 255.0);
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
        public static void FromGrayscale(this double[,] m, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromGrayscale(m, bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
}
