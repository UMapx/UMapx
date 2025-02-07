using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the grayscale filter based on the HSL structure.
    /// <remarks>
    /// The filter discolors the specified part of the image.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class HSLGrayscale : IBitmapFilter
    {
        #region Private data
        private int min;
        private int max;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the grayscale filter based on the HSL structure.
        /// </summary>
        /// <param name="hue">Hue range [0, 359]</param>
        public HSLGrayscale(RangeInt hue)
        {
            Hue = hue;
        }
        /// <summary>
        /// Initializes the grayscale filter based on the HSL structure.
        /// </summary>
        /// <param name="min">Lower bound [0, 359]</param>
        /// <param name="max">Upper bound [0, 359]</param>
        public HSLGrayscale(int min, int max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Gets or sets the hue range [0, 359].
        /// </summary>
        public RangeInt Hue
        {
            get
            {
                return new RangeInt(min, max);
            }
            set
            {
                this.min = value.Min;
                this.max = value.Max;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                HSL hsl; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // This function modifies a given image in order to keep a specific hue
                    // (given too) and to desaturate the rest of the image. This procedure 
                    // originates a image with black and white colormap, excluding the parts
                    // colored with that hue.
                    // Victor Martnez Cagigal, 23/02/2015
                    //
                    // Designed for UMapx.NET by Valery Asiryan, 2018.

                    k = jstride + i * 4;

                    // Convert to hsl:
                    hsl = HSL.FromRGB(p[k + 2], p[k + 1], p[k + 0]);

                    // Getting hue and saturation parameters:
                    float hue = hsl.Hue, saturation = hsl.Saturation;

                    // Applying filter:
                    if (min < max)
                    {
                        hsl.Saturation = (hue > min && hue < max) ? saturation : 0;
                    }
                    else
                    {
                        hsl.Saturation = ((hue > min && hue <= 360) || (hue < max && hue >= 0)) ? saturation : 0;
                    }

                    // Convert to rgb:
                    rgb = hsl.ToRGB;

                    p[k + 0] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                }
            }
            );

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
        }
        #endregion
    }
}
