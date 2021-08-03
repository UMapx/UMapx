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
    public class Resize : IBitmapFilter2
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
            // get source image size
            int width = bmSrc.Width;
            int height = bmSrc.Height;

            int srcStride = bmSrc.Stride;
            int dstOffset = bmData.Stride - 4 * newWidth;
            float xFactor = (float)width / newWidth;
            float yFactor = (float)height / newHeight;

            // do the job
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();

            // coordinates of source points and cooefficiens
            float ox, oy, dx, dy, k1, k2;
            int ox1, oy1, ox2, oy2;
            // destination pixel values
            float r, g, b, a;
            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;
            // temporary pointer
            byte* p;

            // RGB
            for (int y = 0; y < newHeight; y++)
            {
                // Y coordinates
                oy = y * yFactor - 0.5f;
                oy1 = (int)oy;
                dy = oy - oy1;

                for (int x = 0; x < newWidth; x++, dst += 4)
                {
                    // X coordinates
                    ox = x * xFactor - 0.5f;
                    ox1 = (int)ox;
                    dx = ox - ox1;

                    // initial pixel value
                    r = g = b = a = 0;

                    for (int n = -1; n < 3; n++)
                    {
                        // get Y cooefficient
                        k1 = Kernel.Bicubic(dy - n);

                        oy2 = oy1 + n;
                        if (oy2 < 0)
                            oy2 = 0;
                        if (oy2 > ymax)
                            oy2 = ymax;

                        for (int m = -1; m < 3; m++)
                        {
                            // get X cooefficient
                            k2 = k1 * Kernel.Bicubic(m - dx);

                            ox2 = ox1 + m;
                            if (ox2 < 0)
                                ox2 = 0;
                            if (ox2 > xmax)
                                ox2 = xmax;

                            // get pixel of original image
                            p = src + oy2 * srcStride + ox2 * 4;

                            a += k2 * p[3];
                            r += k2 * p[2];
                            g += k2 * p[1];
                            b += k2 * p[0];
                        }
                    }

                    dst[3] = Maths.Byte(a);
                    dst[2] = Maths.Byte(r);
                    dst[1] = Maths.Byte(g);
                    dst[0] = Maths.Byte(b);
                }
                dst += dstOffset;
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
        #endregion
    }
}
