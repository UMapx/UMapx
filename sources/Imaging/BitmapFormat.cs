using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;

namespace UMapx.Imaging
{
    /// <summary>
    /// Used to work with bitmap formats.
    /// </summary>
    public static partial class BitmapFormat
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
            stream.Position = 0;

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
            stream.Position = 0;

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
            stream.Position = 0;

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
            stream.Position = 0;

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
            stream.Position = 0;

            return new Bitmap(stream);
        }
        /// <summary>
        /// Converts BitmapData to Bitmap.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap ToBitmap(this BitmapData bmData)
        {
            int w = bmData.Width;
            int h = bmData.Height;
            var format = bmData.PixelFormat;

            var dst = new Bitmap(w, h, format);
            var rect = new Rectangle(0, 0, w, h);
            var dstData = dst.LockBits(rect, ImageLockMode.WriteOnly, format);

            try
            {
                int bpp = Image.GetPixelFormatSize(format) / 8;
                int rowBytes = checked(w * bpp);

                byte* sBase = (byte*)bmData.Scan0;
                int sStride = bmData.Stride;

                byte* dBase = (byte*)dstData.Scan0;
                int dStride = dstData.Stride;

                for (int y = 0; y < h; y++)
                {
                    byte* s = sBase + (sStride > 0 ? y * sStride : (h - 1 - y) * (-sStride));
                    byte* d = dBase + y * dStride;
                    Buffer.MemoryCopy(s, d, rowBytes, rowBytes);
                }

                return dst;
            }
            finally
            {
                dst.UnlockBits(dstData);
            }
        }
        /// <summary>
        /// Converts Bitmap to 32bpp ARGB format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap To32bpp(this Bitmap b)
        {
            return b.Clone(new Rectangle(0, 0, b.Width, b.Height), PixelFormat.Format32bppArgb);
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
        /// Unblocks Bitmap in system memory.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="bmData">Bitmap data</param>
        public static void Unlock(this Bitmap b, BitmapData bmData)
        {
            b.UnlockBits(bmData);
        }
        #endregion
    }
}
