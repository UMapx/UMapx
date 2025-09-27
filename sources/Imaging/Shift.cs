using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the shift filter.
    /// </summary>
    [Serializable]
    public class Shift : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int x;
        private int y;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the shift filter.
        /// </summary>
        /// <param name="x">Offset value of axis X</param>
        /// <param name="y">Offset value of axis Y</param>
        public Shift(int x, int y)
        {
            X = x;
            Y = y;
        }
        /// <summary>
        /// Initializes the shift filter.
        /// </summary>
        /// <param name="point">A pair of integers representing an ordered pair of X and Y coordinates</param>
        public Shift(PointInt point)
        {
            X = point.X;
            Y = point.Y;
        }
        /// <summary>
        /// Initializes the shift filter.
        /// </summary>
        public Shift() { }
        /// <summary>
        /// Gets or sets the offset value of axis X.
        /// </summary>
        public int X
        {
            get
            {
                return this.x;
            }
            set
            {
                this.x = value;
            }
        }
        /// <summary>
        /// Gets or sets the offset value of axis Y.
        /// </summary>
        public int Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            if (bmData.PixelFormat != PixelFormat.Format32bppArgb || bmSrc.PixelFormat != PixelFormat.Format32bppArgb)
                throw new NotSupportedException("Only support Format32bppArgb pixelFormat");

            int width = bmSrc.Width;
            int height = bmSrc.Height;
            if (bmData.Width != width || bmData.Height != height)
                throw new ArgumentException("Input image sizes must be the same");

            int dx = ((X % width) + width) % width;
            int dy = ((Y % height) + height) % height;
            if (dx == 0 && dy == 0) return;

            int srcStride = bmSrc.Stride;
            int dstStride = bmData.Stride;

            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            byte* pDst = (byte*)bmData.Scan0.ToPointer();

            int rowBytes = checked(width * 4);
            if (rowBytes > srcStride || rowBytes > dstStride)
                throw new ArgumentException("Row size exceeds stride");

            byte* tempRow = stackalloc byte[rowBytes];

            for (int j = 0; j < height; j++)
            {
                int sy = j - dy; if (sy < 0) sy += height;

                byte* srcRow = pSrc + (long)srcStride * sy;
                byte* dstRow = pDst + (long)dstStride * j;

                if (dx == 0)
                {
                    Buffer.MemoryCopy(srcRow, dstRow, dstStride, rowBytes);
                }
                else
                {
                    int splitBytes = dx * 4;
                    int rightBytes = rowBytes - splitBytes;

                    Buffer.MemoryCopy(srcRow + splitBytes, tempRow, rowBytes, rightBytes);
                    Buffer.MemoryCopy(srcRow, tempRow + rightBytes, rowBytes - rightBytes, splitBytes);
                    Buffer.MemoryCopy(tempRow, dstRow, dstStride, rowBytes);
                }
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
