using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the rotation filter.
    /// </summary>
    [Serializable]
    public class Rotate : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private float angle;
        private Color color;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the rotation filter.
        /// </summary>
        /// <param name="angle">Angle</param>
        public Rotate(float angle)
        {
            Angle = angle;
            Color = Color.Transparent;
        }
        /// <summary>
        /// Initializes the rotation filter.
        /// </summary>
        /// <param name="angle">Angle</param>
        /// <param name="color">Background color</param>
        public Rotate(float angle, Color color)
        {
            Angle = angle; Color = color;
        }
        /// <summary>
        /// Gets or sets angle value.
        /// </summary>
        public float Angle
        {
            get
            {
                return this.angle;
            }
            set
            {
                this.angle = value;
            }
        }
        /// <summary>
        /// Gets or sets background color.
        /// </summary>
        public Color Color
        {
            get
            {
                return this.color;
            }
            set
            {
                this.color = value;
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
            float oldXradius = (float)(width - 1) / 2;
            float oldYradius = (float)(height - 1) / 2;

            // get destination image size
            int newWidth = bmData.Width;
            int newHeight = bmData.Height;
            float newXradius = (float)(newWidth - 1) / 2;
            float newYradius = (float)(newHeight - 1) / 2;

            // angle's sine and cosine
            float angleRad = -angle * Maths.Pi / 180.0f;
            float angleCos = Maths.Cos(angleRad);
            float angleSin = Maths.Sin(angleRad);

            int srcStride = bmSrc.Stride;
            int dstOffset = bmData.Stride - newWidth * 4;

            // fill values
            byte fillA = color.A;
            byte fillR = color.R;
            byte fillG = color.G;
            byte fillB = color.B;

            // do the job
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();

            // destination pixel's coordinate relative to image center
            float cx, cy;
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
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                cx = -newXradius;
                for (int x = 0; x < newWidth; x++, dst += 4)
                {
                    // coordinates of source point
                    ox = angleCos * cx + angleSin * cy + oldXradius;
                    oy = -angleSin * cx + angleCos * cy + oldYradius;

                    ox1 = (int)ox;
                    oy1 = (int)oy;

                    // validate source pixel's coordinates
                    if ((ox1 < 0) || (oy1 < 0) || (ox1 >= width) || (oy1 >= height))
                    {
                        // fill destination image with filler
                        dst[3] = fillA;
                        dst[2] = fillR;
                        dst[1] = fillG;
                        dst[0] = fillB;
                    }
                    else
                    {
                        dx = ox - ox1;
                        dy = oy - oy1;

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
                    cx++;
                }
                cy++;
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
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            Apply(bmData);
            BitmapFormat.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
}
