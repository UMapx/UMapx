using System;
using SkiaDrawing;
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
        /// <param name="interpolationMode">Interpolation mode</param>
        public Rotate(float angle, UMapx.Core.InterpolationMode interpolationMode = UMapx.Core.InterpolationMode.Bicubic)
        {
            Angle = angle;
            Color = Color.Transparent;
            InterpolationMode = interpolationMode;
        }
        /// <summary>
        /// Initializes the rotation filter.
        /// </summary>
        /// <param name="angle">Angle</param>
        /// <param name="color">Background color</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        public Rotate(float angle, Color color, UMapx.Core.InterpolationMode interpolationMode = UMapx.Core.InterpolationMode.Bicubic)
        {
            Angle = angle;
            Color = color;
            InterpolationMode = interpolationMode;
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
        /// Gets or sets interpolation mode.
        /// </summary>
        public UMapx.Core.InterpolationMode InterpolationMode { get; set; }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            if (InterpolationMode == UMapx.Core.InterpolationMode.Bicubic)
            {
                ApplyBicubic(bmData, bmSrc);
            }
            else if (InterpolationMode == UMapx.Core.InterpolationMode.Bilinear)
            {
                ApplyBilinear(bmData, bmSrc);
            }
            else if (InterpolationMode == UMapx.Core.InterpolationMode.NearestNeighbor)
            {
                ApplyNearestNeighbor(bmData, bmSrc);
            }
            else
            {
                throw new NotSupportedException();
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
            var Src = (Bitmap)Data.Clone();
            Apply(Data, Src);
            Src.Dispose();
            return;
        }
        #endregion

        #region Private methods

        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyNearestNeighbor(BitmapData bmData, BitmapData bmSrc)
        {
            // get source image size
            int width = bmSrc.Width;
            int height = bmSrc.Height;
            double oldXradius = (double)(width - 1) / 2;
            double oldYradius = (double)(height - 1) / 2;

            // get destination image size
            int newWidth = bmData.Width;
            int newHeight = bmData.Height;
            double newXradius = (double)(newWidth - 1) / 2;
            double newYradius = (double)(newHeight - 1) / 2;

            // angle's sine and cosine
            double angleRad = -angle * Math.PI / 180;
            double angleCos = Math.Cos(angleRad);
            double angleSin = Math.Sin(angleRad);

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
            double cx, cy;
            // source pixel's coordinates
            int ox, oy;
            // temporary pointer
            byte* p;

            // check pixel format
            // ARGB
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                cx = -newXradius;
                for (int x = 0; x < newWidth; x++, dst += 4)
                {
                    // coordinate of the nearest point
                    ox = (int)(angleCos * cx + angleSin * cy + oldXradius);
                    oy = (int)(-angleSin * cx + angleCos * cy + oldYradius);

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
                        p = src + oy * srcStride + ox * 4;

                        dst[3] = p[3];
                        dst[2] = p[2];
                        dst[1] = p[1];
                        dst[0] = p[0];
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
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyBilinear(BitmapData bmData, BitmapData bmSrc)
        {
            // get source image size
            int width = bmSrc.Width;
            int height = bmSrc.Height;
            double oldXradius = (double)(width - 1) / 2;
            double oldYradius = (double)(height - 1) / 2;

            // get destination image size
            int newWidth = bmData.Width;
            int newHeight = bmData.Height;
            double newXradius = (double)(newWidth - 1) / 2;
            double newYradius = (double)(newHeight - 1) / 2;

            // angle's sine and cosine
            double angleRad = -angle * Math.PI / 180;
            double angleCos = Math.Cos(angleRad);
            double angleSin = Math.Sin(angleRad);

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
            double cx, cy;
            // coordinates of source points
            double ox, oy, tx, ty, dx1, dy1, dx2, dy2;
            int ox1, oy1, ox2, oy2;
            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;
            // temporary pointers
            byte* p1, p2, p3, p4;

            // RGB
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                // do some pre-calculations of source points' coordinates
                // (calculate the part which depends on y-loop, but does not
                // depend on x-loop)
                tx = angleSin * cy + oldXradius;
                ty = angleCos * cy + oldYradius;

                cx = -newXradius;
                for (int x = 0; x < newWidth; x++, dst += 4)
                {
                    // coordinates of source point
                    ox = tx + angleCos * cx;
                    oy = ty - angleSin * cx;

                    // top-left coordinate
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
                        // bottom-right coordinate
                        ox2 = (ox1 == xmax) ? ox1 : ox1 + 1;
                        oy2 = (oy1 == ymax) ? oy1 : oy1 + 1;

                        if ((dx1 = ox - (float)ox1) < 0)
                            dx1 = 0;
                        dx2 = 1.0f - dx1;

                        if ((dy1 = oy - (float)oy1) < 0)
                            dy1 = 0;
                        dy2 = 1.0f - dy1;

                        // get four points
                        p1 = p2 = src + oy1 * srcStride;
                        p1 += ox1 * 4;
                        p2 += ox2 * 4;

                        p3 = p4 = src + oy2 * srcStride;
                        p3 += ox1 * 4;
                        p4 += ox2 * 4;

                        // interpolate using 4 points
                        // alpha
                        dst[3] = (byte)(
                            dy2 * (dx2 * p1[3] + dx1 * p2[3]) +
                            dy1 * (dx2 * p3[3] + dx1 * p4[3]));

                        // red
                        dst[2] = (byte)(
                            dy2 * (dx2 * p1[2] + dx1 * p2[2]) +
                            dy1 * (dx2 * p3[2] + dx1 * p4[2]));

                        // green
                        dst[1] = (byte)(
                            dy2 * (dx2 * p1[1] + dx1 * p2[1]) +
                            dy1 * (dx2 * p3[1] + dx1 * p4[1]));

                        // blue
                        dst[0] = (byte)(
                            dy2 * (dx2 * p1[0] + dx1 * p2[0]) +
                            dy1 * (dx2 * p3[0] + dx1 * p4[0]));
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
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyBicubic(BitmapData bmData, BitmapData bmSrc)
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

        #endregion
    }
}
