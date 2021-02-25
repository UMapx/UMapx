using System;
using System.Drawing;
using System.Drawing.Imaging;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the rotation filter.
    /// </summary>
    [Serializable]
    public class Rotate : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private double angle;
        private Color color;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the rotation filter.
        /// </summary>
        /// <param name="angle">Angle</param>
        /// <param name="color">Background color</param>
        public Rotate(double angle, Color color)
        {
            Angle = angle; Color = color;
        }
        /// <summary>
        /// Gets or sets angle value.
        /// </summary>
        public double Angle
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
            // exception!
            if (bmSrc.Width != bmData.Width || bmSrc.Height != bmData.Height)
                throw new Exception("Input image sizes must be the same");

            // get source image size
            int width = bmSrc.Width, height = bmSrc.Height, stride = bmSrc.Stride;
            double xradius = (width - 1) / 2.0, yradius = (height - 1) / 2.0;

            // angle's sine and cosine
            double angleRad = -angle * Math.PI / 180;
            double angleCos = Math.Cos(angleRad);
            double angleSin = Math.Sin(angleRad);

            // destination pixel's coordinate relative to image center
            double cx, cy = -yradius;

            // fill values
            byte fillA = color.A;
            byte fillR = color.R;
            byte fillG = color.G;
            byte fillB = color.B;

            // do the job
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* p;

            // source pixel's coordinates
            int ox, oy, y, x;

            for (y = 0; y < height; y++)
            {
                cx = -xradius;
                for (x = 0; x < width; x++, dst += 4)
                {
                    // coordinate of the nearest point
                    ox = (int)(angleCos * cx + angleSin * cy + xradius);
                    oy = (int)(-angleSin * cx + angleCos * cy + yradius);

                    // validate source pixel's coordinates
                    if ((ox < 0) || (oy < 0) || (ox >= width) || (oy >= height))
                    {
                        // fill destination image with filler
                        dst[3] = fillA;
                        dst[2] = fillR;
                        dst[1] = fillG;
                        dst[0] = fillB;
                    }
                    else
                    {
                        // fill destination image with pixel from source image
                        p = src + oy * stride + ox * 4;

                        dst[3] = p[3];
                        dst[2] = p[2];
                        dst[1] = p[1];
                        dst[0] = p[0];
                    }
                    cx++;
                }
                cy++;
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
