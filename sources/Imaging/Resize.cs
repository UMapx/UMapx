using System;
using SkiaDrawing;
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
        /// <param name="interpolationMode">Interpolation mode</param>
        public Resize(int width = 512, int height = 512, UMapx.Core.InterpolationMode interpolationMode = UMapx.Core.InterpolationMode.Bicubic)
        {
            Size = new SizeInt(width, height);
            InterpolationMode = interpolationMode;
        }
        /// <summary>
        /// Initializes the resize filter.
        /// </summary>
        /// <param name="size">Size</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        public Resize(SizeInt size, UMapx.Core.InterpolationMode interpolationMode = UMapx.Core.InterpolationMode.Bicubic)
        {
            Size = size;
            InterpolationMode = interpolationMode;
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

            int srcStride = bmSrc.Stride;
            int dstStride = bmData.Stride;
            double xFactor = (double)width / newWidth;
            double yFactor = (double)height / newHeight;

            // do the job
            byte* baseSrc = (byte*)bmSrc.Scan0.ToPointer();
            byte* baseDst = (byte*)bmData.Scan0.ToPointer();

            // for each line
            for (int y = 0; y < newHeight; y++)
            {
                byte* dst = baseDst + dstStride * y;
                byte* src = baseSrc + srcStride * ((int)(y * yFactor));
                byte* p;

                // for each pixel
                for (int x = 0; x < newWidth; x++)
                {
                    p = src + 4 * ((int)(x * xFactor));

                    for (int i = 0; i < 4; i++, dst++, p++)
                    {
                        *dst = *p;
                    }
                }
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

            int srcStride = bmSrc.Stride;
            int dstOffset = bmData.Stride - 4 * newWidth;
            double xFactor = (double)width / newWidth;
            double yFactor = (double)height / newHeight;

            // do the job
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();

            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;
            // temporary pointers

            // for each line
            for (int y = 0; y < newHeight; y++)
            {
                // Y coordinates
                double oy = (double)y * yFactor;
                int oy1 = (int)oy;
                int oy2 = (oy1 == ymax) ? oy1 : oy1 + 1;
                double dy1 = oy - (double)oy1;
                double dy2 = 1.0 - dy1;

                // get temp pointers
                byte* tp1 = src + oy1 * srcStride;
                byte* tp2 = src + oy2 * srcStride;

                // for each pixel
                for (int x = 0; x < newWidth; x++)
                {
                    // X coordinates
                    double ox = (double)x * xFactor;
                    int ox1 = (int)ox;
                    int ox2 = (ox1 == xmax) ? ox1 : ox1 + 1;
                    double dx1 = ox - (double)ox1;
                    double dx2 = 1.0 - dx1;

                    // get four points
                    byte* p1 = tp1 + ox1 * 4;
                    byte* p2 = tp1 + ox2 * 4;
                    byte* p3 = tp2 + ox1 * 4;
                    byte* p4 = tp2 + ox2 * 4;

                    // interpolate using 4 points
                    for (int i = 0; i < 4; i++, dst++, p1++, p2++, p3++, p4++)
                    {
                        *dst = (byte)(
                            dy2 * (dx2 * (*p1) + dx1 * (*p2)) +
                            dy1 * (dx2 * (*p3) + dx1 * (*p4)));
                    }
                }
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


        #endregion
    }
}
