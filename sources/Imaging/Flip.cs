using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the flip filter.
    /// </summary>
    [Serializable]
    public class Flip : IBitmapFilter
    {
        #region Private data
        private bool x = true;
        private bool y = true;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the flip filter.
        /// </summary>
        public Flip() { }
        /// <summary>
        /// Initializes the flip filter.
        /// </summary>
        /// <param name="x">Flip X</param>
        /// <param name="y">Flip Y</param>
        public Flip(bool x, bool y)
        {
            X = x;
            Y = y;
        }
        /// <summary>
        /// Gets or sets flip X.
        /// </summary>
        public bool X
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
        /// Gets or sets flip Y.
        /// </summary>
        public bool Y
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
        public unsafe void Apply(BitmapData bmData)
        {
            #region Data
            int pixel = 4, start = 0;
            int width = bmData.Width, height = bmData.Height;
            int offset = width * pixel, stride = bmData.Stride;
            #endregion

            #region FlipY
            // Flip X:
            if (this.x == true)
            {
                byte* p = (byte*)bmData.Scan0.ToPointer();
                byte* pSrc = (byte*)bmData.Scan0.ToPointer();
                int s0 = stride - (width >> 1) * pixel;
                int s1 = stride + (width >> 1) * pixel;
                pSrc += (width - 1) * pixel;
                int l, w2;
                byte b;

                for (int k = 0; k < height; k++)
                {
                    l = 0; w2 = (width >> 1);

                    while (l < w2)
                    {
                        b = p[2];
                        p[2] = pSrc[2];
                        pSrc[2] = b;

                        b = p[1];
                        p[1] = pSrc[1];
                        pSrc[1] = b;

                        b = *p;
                        *p = *pSrc;
                        *pSrc = b;

                        l++;
                        p += 4;
                        pSrc -= 4;
                    }
                    p += s0;
                    pSrc += s1;
                }
            }
            #endregion

            #region FlipX
            // Flip Y:
            if (this.y == true)
            {
                byte* p = (byte*)bmData.Scan0.ToPointer();
                byte* pSrc = (byte*)bmData.Scan0.ToPointer();
                pSrc += (height - 1) * stride;
                int offset2 = stride - width * pixel;
                int m = 0, n = start;
                int h2 = (height >> 1);
                byte b2;

                while (m < h2)
                {
                    n = start;
                    while (n < offset)
                    {
                        b2 = *p;
                        *p = *pSrc;
                        *pSrc = b2;
                        n++;
                        p++;
                        pSrc++;
                    }
                    p += offset2;
                    pSrc += offset2 - stride - stride;
                    m++;
                }
            }
            #endregion
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
        }
        #endregion
    }
}
