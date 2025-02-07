using System;
using SkiaDrawing;
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
            // get source image:
            int width = bmSrc.Width;
            int height = bmSrc.Height;

            // images strides:
            int srcStride = bmSrc.Stride;
            int dstStride = bmData.Stride;

            // do the job:
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            byte* pDst = (byte*)bmData.Scan0.ToPointer();
            byte* dst, src, p;

            // source pixel's coordinates from rectangle:
            int startX = Maths.Range(rectangle.X, 0, width);
            int startY = Maths.Range(rectangle.Y, 0, height);
            int endX = Maths.Range(rectangle.Width - startX, 0, width);
            int endY = Maths.Range(rectangle.Height - startY, 0, height);

            // pixel offsets:
            int x, y, i;

            // for each line:
            for (y = 0; y < endY; y++)
            {
                dst = pDst + dstStride * y;
                src = pSrc + srcStride * (y + startX);

                // for each pixel:
                for (x = 0; x < endX; x++)
                {
                    p = src + 4 * (x + startX);

                    for (i = 0; i < 4; i++, dst++, p++)
                    {
                        *dst = *p;
                    }
                }
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
    }
}
