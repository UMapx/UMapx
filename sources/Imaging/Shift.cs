using System;
using SkiaDrawing;
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
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // image properties:
            int width = bmSrc.Width;
            int height = bmSrc.Height;
            int stride = bmSrc.Stride;

            // exception!
            if (bmData.Width != width ||
                bmData.Height != height)
                throw new Exception("Input image sizes must be the same");

            // applying only for X:
            if (x != 0 && y == 0)
                ShiftX(bmData, bmSrc, width, height, stride);

            // applying only for Y:
            else if (y != 0 && x == 0)
                ShiftY(bmData, bmSrc, width, height, stride);

            // applying for two sides:
            else
            {
                ShiftX(bmSrc, bmData, width, height, stride);
                ShiftY(bmData, bmSrc, width, height, stride);
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

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="stride">Stride</param>
        private unsafe void ShiftY(BitmapData bmData, BitmapData bmSrc, int width, int height, int stride)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)(void*)bmSrc.Scan0.ToPointer();
            int sX = 0, sY = 0, front = 0, reverse = 0, f = 1;
            int i, j;

            for (j = 0; j < height; j++)
            {
                sY = j + y; // Offset of axis Y

                for (i = 0; i < width; i++, p += 4)
                {
                    sX = i; // Offset of axis X

                    if (sY < height && sY >= 0)
                    {
                        front = Math.Abs(sX * 4 + sY * stride);
                        p[0] = pSrc[front];
                        p[1] = pSrc[front + 1];
                        p[2] = pSrc[front + 2];
                        p[3] = pSrc[front + 3];
                    }
                    else
                    {
                        f = (sY < 0) ? -1 : 1;

                        reverse = Math.Abs(sX * 4 + (sY - f * height) * stride);
                        p[0] = pSrc[reverse];
                        p[1] = pSrc[reverse + 1];
                        p[2] = pSrc[reverse + 2];
                        p[3] = pSrc[reverse + 3];
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="stride">Stride</param>
        private unsafe void ShiftX(BitmapData bmData, BitmapData bmSrc, int width, int height, int stride)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)(void*)bmSrc.Scan0.ToPointer();
            int sX = 0, sY = 0, front = 0, reverse = 0, f = 1;
            int i, j;

            for (j = 0; j < height; j++)
            {
                sY = j; // Offset of axis Y

                for (i = 0; i < width; i++, p += 4)
                {
                    sX = i + x; // Offset of axis X

                    if (sX < width && sX >= 0)
                    {
                        front = Math.Abs(sX * 4 + sY * stride);
                        p[0] = pSrc[front];
                        p[1] = pSrc[front + 1];
                        p[2] = pSrc[front + 2];
                        p[3] = pSrc[front + 3];
                    }
                    else
                    {
                        f = (sX < 0) ? -1 : 1;

                        reverse = Math.Abs((sX - f * width) * 4 + sY * stride);
                        p[0] = pSrc[reverse];
                        p[1] = pSrc[reverse + 1];
                        p[2] = pSrc[reverse + 2];
                        p[3] = pSrc[reverse + 3];
                    }
                }
            }
            return;
        }
        #endregion
    }
}
