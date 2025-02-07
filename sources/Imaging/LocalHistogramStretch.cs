using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the local histogram stretch filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.academia.edu/7629047/Image_enhancement_by_local_histogram_stretching
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LocalHistogramStretch : IBitmapFilter
    {
        #region Private data
        private readonly BoxBlur gb = new BoxBlur();
        private readonly Erosion er = new Erosion();
        private readonly Dilatation di = new Dilatation();
        private Space space;
        private float contrast;
        private bool smoothing;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local histogram stretch filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [0, 1]</param>
        /// <param name="smoothing">Smoothing</param>
        public LocalHistogramStretch(int radius, Space space, float contrast = 0.5f, bool smoothing = true)
        {
            Size = new SizeInt(radius, radius);
            Space = space;
            Smoothing = smoothing;
            Contrast = contrast;
        }
        /// <summary>
        /// Initializes the local histogram stretch filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [0, 1]</param>
        /// <param name="smoothing">Smoothing</param>
        public LocalHistogramStretch(int width, int height, Space space, float contrast = 0.5f, bool smoothing = true)
        {
            Size = new SizeInt(width, height);
            Space = space;
            Smoothing = smoothing;
            Contrast = contrast;
        }
        /// <summary>
        /// Initializes the local histogram stretch filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [0, 1]</param>
        /// <param name="smoothing">Smoothing</param>
        public LocalHistogramStretch(SizeInt size, Space space, float contrast = 0.5f, bool smoothing = true)
        {
            Size = size;
            Space = space;
            Smoothing = smoothing;
            Contrast = contrast;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return er.Size;
            }
            set
            {
                er.Size = value; // local min  filter -> r,
                di.Size = value; // local max  filter -> r,
                gb.Size =        // local mean filter -> 2r.
                    new SizeInt( // we got only even values for radius.
                        value.Width * 2,
                        value.Height * 2);
            }
        }
        /// <summary>
        /// Gets or sets the color space.
        /// </summary>
        public Space Space
        {
            get
            {
                return this.space;
            }
            set
            {
                this.space = value;
            }
        }
        /// <summary>
        /// Gets or sets the contrast value [0, 1].
        /// </summary>
        public float Contrast
        {
            get
            {
                return this.contrast;
            }
            set
            {
                this.contrast = Maths.Float(value);
            }
        }
        /// <summary>
        /// Smoothing.
        /// </summary>
        public bool Smoothing
        {
            get
            {
                return this.smoothing;
            }
            set
            {
                this.smoothing = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            Bitmap current = BitmapFormat.Bitmap(bmData);
            Bitmap Max = (Bitmap)current.Clone();
            Bitmap Min = (Bitmap)Max.Clone();

            di.Apply(Max); er.Apply(Min);

            BitmapData bmMax = BitmapFormat.Lock32bpp(Max);
            BitmapData bmMin = BitmapFormat.Lock32bpp(Min);

            if (smoothing)
                gb.Apply(bmMax); gb.Apply(bmMin);

            Apply(bmData, bmMax, bmMin);

            BitmapFormat.Unlock(Max, bmMax);
            BitmapFormat.Unlock(Min, bmMin);

            Max.Dispose(); Min.Dispose(); current.Dispose();
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Max = (Bitmap)Data.Clone();
            Bitmap Min = (Bitmap)Data.Clone();

            di.Apply(Max); er.Apply(Min);

            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            BitmapData bmMax = BitmapFormat.Lock32bpp(Max);
            BitmapData bmMin = BitmapFormat.Lock32bpp(Min);

            if (smoothing)
                gb.Apply(bmMax); gb.Apply(bmMin);

            Apply(bmData, bmMax, bmMin);

            BitmapFormat.Unlock(Data, bmData);
            BitmapFormat.Unlock(Max, bmMax);
            BitmapFormat.Unlock(Min, bmMin);

            Max.Dispose(); Min.Dispose();
            return;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void Apply(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            // filter
            switch (space)
            {
                case Imaging.Space.HSB:
                    ApplyHSB(bmData, bmMax, bmMin);
                    break;
                case Imaging.Space.HSL:
                    ApplyHSL(bmData, bmMax, bmMin);
                    break;
                case Imaging.Space.YCbCr:
                    ApplyYCbCr(bmData, bmMax, bmMin);
                    break;
                case Imaging.Space.RGB:
                    ApplyRGB(bmData, bmMax, bmMin);
                    break;
                default:
                    ApplyGrayscale(bmData, bmMax, bmMin);
                    break;
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyRGB(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            float required = 1.0f - this.contrast;

            Parallel.For(0, height, j =>
            {
                int i, k, k1, k2, q, v, jstride = j * stride;
                float mag, max, min;
                float num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    // Local function:
                    for (q = 0; q < 3; q++)
                    {
                        v = k + q;

                        mag = p[v] / 255.0f;
                        max = pMax[v] / 255.0f;
                        min = pMin[v] / 255.0f;

                        num1 = max - min;

                        if (num1 < required)
                        {
                            num2 = min + (required - num1) * min / (num1 - 1f);
                            min = Maths.Float(num2);
                            max = Maths.Float(num2 + required);
                        }

                        num1 = max - min;
                        num3 = mag - min;

                        if (num1 > 0)
                        {
                            p[v] = Maths.Byte(255 * num3 / num1);
                        }
                    }
                    // end local function.
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyHSB(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            float required = 1.0f - this.contrast;

            Parallel.For(0, height, j =>
            {
                HSB imag; HSB imax; HSB imin; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;
                float mag, max, min;
                float num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    imag = HSB.FromRGB(p[k2], p[k1], p[k]);
                    imax = HSB.FromRGB(pMax[k2], pMax[k1], pMax[k]);
                    imin = HSB.FromRGB(pMin[k2], pMin[k1], pMin[k]);

                    mag = imag.Brightness;
                    max = imax.Brightness;
                    min = imin.Brightness;

                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Float(num2);
                        max = Maths.Float(num2 + required);
                    }

                    num1 = max - min;
                    num3 = mag - min;

                    if (num1 > 0)
                    {
                        imag.Brightness = num3 / num1;
                        rgb = imag.ToRGB;
                        p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                    }
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyHSL(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            float required = 1.0f - this.contrast;

            Parallel.For(0, height, j =>
            {
                HSL imag; HSL imax; HSL imin; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;
                float mag, max, min;
                float num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    imag = HSL.FromRGB(p[k2], p[k1], p[k]);
                    imax = HSL.FromRGB(pMax[k2], pMax[k1], pMax[k]);
                    imin = HSL.FromRGB(pMin[k2], pMin[k1], pMin[k]);

                    mag = imag.Lightness;
                    max = imax.Lightness;
                    min = imin.Lightness;

                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Float(num2);
                        max = Maths.Float(num2 + required);
                    }

                    num1 = max - min;
                    num3 = mag - min;

                    if (num1 > 0)
                    {
                        imag.Lightness = num3 / num1;
                        rgb = imag.ToRGB;
                        p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                    }
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyYCbCr(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            float required = 1.0f - this.contrast;

            Parallel.For(0, height, j =>
            {
                YCbCr imag; YCbCr imax; YCbCr imin; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;
                float mag, max, min;
                float num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    imag = YCbCr.FromRGB(p[k2], p[k1], p[k]);
                    imax = YCbCr.FromRGB(pMax[k2], pMax[k1], pMax[k]);
                    imin = YCbCr.FromRGB(pMin[k2], pMin[k1], pMin[k]);

                    mag = imag.Y;
                    max = imax.Y;
                    min = imin.Y;

                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Float(num2);
                        max = Maths.Float(num2 + required);
                    }

                    num1 = max - min;
                    num3 = mag - min;

                    if (num1 > 0)
                    {
                        imag.Y = num3 / num1;
                        rgb = imag.ToRGB;
                        p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                    }
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyGrayscale(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            float required = 1.0f - this.contrast;

            Parallel.For(0, height, j =>
            {
                int i, k, k1, k2, v, jstride = j * stride;
                float mag, max, min;
                float num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    // Local function:
                    v = k;

                    mag = RGB.Average(p[v], p[v + 1], p[v + 2]) / 255.0f;
                    max = RGB.Average(pMax[v], pMax[v + 1], pMax[v + 2]) / 255.0f;
                    min = RGB.Average(pMin[v], pMin[v + 1], pMin[v + 2]) / 255.0f;

                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Float(num2);
                        max = Maths.Float(num2 + required);
                    }

                    num1 = max - min;
                    num3 = mag - min;

                    if (num1 > 0)
                    {
                        p[v] = p[v + 1] = p[v + 2] = Maths.Byte(255 * num3 / num1);
                    }
                    // end local function.
                }
            }
            );
        }
        #endregion
    }
}
