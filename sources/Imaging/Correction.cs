using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the mask correction filter.
    /// </summary>
    [Serializable]
    public class Correction : Rebuilder, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Values.
        /// </summary>
        protected float[] values;
        /// <summary>
        /// Color space.
        /// </summary>
        protected Space space;
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the mask correction filter.
        /// </summary>
        /// <param name="values">Mask array</param>
        /// <param name="space">Color space</param>
        public Correction(float[] values, Space space)
        {
            Values = values;
            Space = space;
        }
        /// <summary>
        /// Initializes the mask correction filter.
        /// </summary>
        public Correction()
        {
            Values = new float[256];
            Space = Imaging.Space.RGB;
        }
        /// <summary>
        /// Gets or sets the mask array.
        /// </summary>
        public float[] Values
        {
            get
            {
                return this.values;
            }
            set
            {
                if (value.Length != 256)
                    throw new Exception("Mask size should be  256");

                this.values = value;
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
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            // rebuild?
            if (rebuild == true)
            {
                this.Rebuild(); this.rebuild = false;
            }

            // filter
            switch (space)
            {
                case Imaging.Space.HSB:
                    ApplyHSB(bmData);
                    break;
                case Imaging.Space.HSL:
                    ApplyHSL(bmData);
                    break;
                case Imaging.Space.YCbCr:
                    ApplyYCbCr(bmData);
                    break;
                case Imaging.Space.RGB:
                    ApplyRGB(bmData);
                    break;
                default:
                    ApplyGrayscale(bmData);
                    break;
            }
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

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyRGB(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                RGB rgb;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    rgb = new RGB(p[k + 2], p[k + 1], p[k]);

                    rgb.Red = Maths.Byte(values[rgb.Red] * length);
                    rgb.Green = Maths.Byte(values[rgb.Green] * length);
                    rgb.Blue = Maths.Byte(values[rgb.Blue] * length);

                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyHSL(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                HSL hsl; RGB rgb;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    hsl = HSL.FromRGB(p[k + 2], p[k + 1], p[k]);
                    hsl.Lightness = values[(int)(hsl.Lightness * length)];
                    rgb = hsl.ToRGB;
                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyHSB(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                HSB hsb; RGB rgb;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    hsb = HSB.FromRGB(p[k + 2], p[k + 1], p[k]);
                    hsb.Brightness = values[(int)(hsb.Brightness * length)];
                    rgb = hsb.ToRGB;
                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyYCbCr(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                YCbCr ycbcr; RGB rgb;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    ycbcr = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k]);
                    ycbcr.Y = values[(int)(ycbcr.Y * length)];
                    rgb = ycbcr.ToRGB;
                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyGrayscale(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                int luma;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    luma = RGB.Average(p[k + 2], p[k + 1], p[k]);
                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(values[luma] * length);
                }
            }
            );

            return;
        }
        #endregion
    }
}
