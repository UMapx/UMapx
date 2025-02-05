//using System.Drawing;
//using System.Drawing.Imaging;

using System;
using System.IO;
using SkiaDrawing;


namespace UMapx.Imaging
{
    /// <summary>
    /// Uses to work with bitmap formats.
    /// </summary>
    public static class BitmapFormat
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
        public static BitmapData Lock24bpp(Bitmap b)
        {
            return b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
        }
        /// <summary>
        /// Blocks Bitmap in system memory.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap data</returns>
        public static BitmapData Lock8bpp(Bitmap b)
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
        }
        #endregion

        #region Bitmap voids
        /// <summary>
        /// Returns new bitmap.
        /// </summary>
        /// <param name="size">Size</param>
        /// <param name="color">Color</param>
        /// <returns>Bitmap</returns>
        public static Bitmap CreateBitmap(this Size size, Color color)
        {
            var bitmap = new Bitmap(size.Width, size.Height);
            using var graphics = Graphics.FromImage(bitmap);
            graphics.Clear(color);
            return bitmap;
        }
        #endregion
    }
}
