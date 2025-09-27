using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Used to work with depth matrices.
    /// </summary>
    public static partial class DepthMatrix
    {
        #region Depth convert components
        /// <summary>
        /// Converts Bitmap into ushort matrix.
        /// </summary>
        /// <param name="bitmap">Bitmap</param>
        /// <returns>Depth</returns>
        /// <remarks>It locks bitmap in 24bpp RGB format.</remarks>
        public unsafe static ushort[,] ToDepth(this Bitmap bitmap)
        {
            var width = bitmap.Width;
            var height = bitmap.Height;
            var rectangle = new Rectangle(0, 0, width, height);
            var bmData = bitmap.LockBits(rectangle, ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            var output = ToDepth(bmData);
            bitmap.Unlock(bmData);
            return output;
        }
        /// <summary>
        /// Converts Bitmap data into ushort matrix.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Depth</returns>
        public unsafe static ushort[,] ToDepth(this BitmapData bmData)
        {
            if (bmData.PixelFormat != PixelFormat.Format24bppRgb && bmData.PixelFormat != PixelFormat.Format32bppArgb)
                throw new NotSupportedException("Depth conversion supports only Format24bppRgb and Format32bppArgb pixel formats.");

            var p = Image.GetPixelFormatSize(bmData.PixelFormat) / 8;
            var width = bmData.Width;
            var height = bmData.Height;
            var stride = bmData.Stride;
            var src = (byte*)bmData.Scan0.ToPointer();
            var output = new ushort[height, width];

            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    var k = x * p + y * stride;
                    output[y, x] = (ushort)RGB.Average(src[k + 2], src[k + 1], src[k + 0]);
                    // ignore alpha channel
                }
            }

            return output;
        }
        /// <summary>
        /// Converts ushort matrix into Bitmap.
        /// </summary>
        /// <param name="depth">Matrix</param>
        /// <returns>Bitmap</returns>
        /// <remarks>It returns bitmap in 24bpp RGB format.</remarks>
        public unsafe static Bitmap FromDepth(this ushort[,] depth)
        {
            var width = depth.GetLength(1);
            var height = depth.GetLength(0);
            var rectangle = new Rectangle(0, 0, width, height);
            var bitmap = new Bitmap(width, height);
            var bmData = bitmap.LockBits(rectangle, ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            var dst = (byte*)bmData.Scan0.ToPointer();
            var stride = bmData.Stride;

            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    var k = x * 3 + y * stride;
                    dst[k + 0] = dst[k + 1] = dst[k + 2] = Maths.Byte((float)depth[y, x] / byte.MaxValue);
                    // ignore alpha channel
                }
            }

            bitmap.Unlock(bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts the depth to the matrix.
        /// </summary>
        /// <param name="depth">Depth</param>
        /// <returns>Matrix</returns>
        public static float[,] ToFloat(this ushort[,] depth)
        {
            int h = depth.GetLength(0);
            int w = depth.GetLength(1);
            float[,] output = new float[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    output[i, j] = depth[i, j] / (float)ushort.MaxValue;
                }
            }

            return output;
        }
        /// <summary>
        /// Converts the ushort matrix to the float matrix.
        /// </summary>
        /// <param name="depth">Matrix</param>
        /// <returns>Matrix</returns>
        public static ushort[,] FromFloat(this float[,] depth)
        {
            int h = depth.GetLength(0);
            int w = depth.GetLength(1);
            ushort[,] output = new ushort[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    output[i, j] = (ushort)(depth[i, j] * ushort.MaxValue);
                }
            }

            return output;
        }
        #endregion
    }
}
