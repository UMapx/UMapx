using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the crop filter.
    /// </summary>
    [Serializable]
    public class Crop : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        Rectangle rectangle;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the crop filter.
        /// </summary>
        /// <param name="rectangle">Rectangle</param>
        public Crop(Rectangle rectangle)
        {
            Rectangle = rectangle;
        }
        /// <summary>
        /// Initializes the crop filter.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        public Crop(int x, int y, int width, int height)
        {
            Rectangle = new Rectangle(x, y, width, height);
        }
        /// <summary>
        /// Gets or sets rectangle.
        /// </summary>
        public Rectangle Rectangle
        {
            get
            {
                return rectangle;
            }
            set
            {
                rectangle = value;
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

            int srcW = bmSrc.Width;
            int srcH = bmSrc.Height;

            int rx = rectangle.X;
            int ry = rectangle.Y;
            int rw = rectangle.Width;
            int rh = rectangle.Height;

            if (rw < 0) { rx += rw; rw = -rw; }
            if (rh < 0) { ry += rh; rh = -rh; }

            int startX = Maths.Range(rx, 0, srcW);
            int startY = Maths.Range(ry, 0, srcH);

            int cropW = Math.Max(0, Math.Min(rw, srcW - startX));
            int cropH = Math.Max(0, Math.Min(rh, srcH - startY));

            if (bmData.Width != cropW || bmData.Height != cropH)
                throw new ArgumentException("Destination BitmapData size must match crop rectangle");

            int srcStride = bmSrc.Stride;
            int dstStride = bmData.Stride;

            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            byte* pDst = (byte*)bmData.Scan0.ToPointer();

            int rowBytes = checked(cropW * 4);

            for (int y = 0; y < cropH; y++)
            {
                byte* srcRow = pSrc + srcStride * (y + startY) + 4 * startX;
                byte* dstRow = pDst + dstStride * y;

                Buffer.MemoryCopy(srcRow, dstRow, dstStride, rowBytes);
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
