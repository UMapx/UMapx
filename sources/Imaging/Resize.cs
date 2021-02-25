using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the resize filter.
    /// </summary>
    [Serializable]
    public class Resize : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        int newWidth;
        int newHeight;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the resize filter.
        /// </summary>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        public Resize(int width = 512, int height = 512)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the resize filter.
        /// </summary>
        /// <param name="size">Size</param>
        public Resize(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Gets or sets image size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(newWidth, newHeight);
            }
            set
            {
                this.newWidth = value.Width;
                this.newHeight = value.Height;
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
            int width = Maths.Range(newWidth, 0, bmData.Width);
            int height = Maths.Range(newHeight, 0, bmData.Height);
            int srcStride = bmSrc.Stride, dstStride = bmData.Stride;
            double xFactor = (double)bmSrc.Width / newWidth;
            double yFactor = (double)bmSrc.Height / newHeight;

            // do the job:
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            byte* pDst = (byte*)bmData.Scan0.ToPointer();
            byte* dst, src, p;

            // source pixel's coordinates
            int x, y, i;

            // for each line
            for (y = 0; y < height; y++)
            {
                dst = pDst + dstStride * y;
                src = pSrc + srcStride * ((int)(y * yFactor));

                // for each pixel
                for (x = 0; x < width; x++)
                {
                    p = src + 4 * ((int)(x * xFactor));

                    for (i = 0; i < 4; i++, dst++, p++)
                    {
                        *dst = *p;
                    }
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
            Bitmap Src = (Bitmap)BitmapFormat.Bitmap(bmData).Clone();
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Src, bmSrc);
            Src.Dispose();
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            Apply(bmData);
            BitmapFormat.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
}
