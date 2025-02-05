using System;
using UMapx.Core;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the linear operation filter.
    /// <remarks>
    /// This filter works UMapxing to the following algorithm: C (x, y) = a * A (x, y) + b * B (x, y), where A, B are the original images,
    /// a, b are the coefficients.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Operation : IBitmapFilter2
    {
        #region Private data
        private float a;
        private float b;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the linear operation filter.
        /// </summary>
        /// <param name="a">First image coefficient</param>
        /// <param name="b">Second image coefficient</param>
        public Operation(float a, float b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the first image coefficient.
        /// </summary>
        public float A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the second image coefficient.
        /// </summary>
        public float B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;

            for (x = 0; x < width; x++)
            {
                for (y = 0; y < height; y++, p += 4, pSrc += 4)
                {
                    p[2] = Maths.Byte(a * p[2] + b * pSrc[2]);
                    p[1] = Maths.Byte(a * p[1] + b * pSrc[1]);
                    p[0] = Maths.Byte(a * p[0] + b * pSrc[0]);
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
        }
        #endregion

        #region Public static methods
        /// <summary>
        /// Addition filter.
        /// </summary>
        public static Operation Addition
        {
            get
            {
                return new Operation(1, 1);
            }
        }
        /// <summary>
        /// Subtraction filter.
        /// </summary>
        public static Operation Subtraction
        {
            get
            {
                return new Operation(1, -1);
            }
        }
        /// <summary>
        /// Averaging filter.
        /// </summary>
        public static Operation Averaging
        {
            get
            {
                return new Operation(0.5f, 0.5f);
            }
        }
        #endregion
    }
}
