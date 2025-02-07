using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the local mask correction filter.
    /// </summary>
    [Serializable]
    public class LocalCorrection : Rebuilder, IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Box blur filter.
        /// </summary>
        protected BoxBlur gb;
        /// <summary>
        /// Contrast.
        /// </summary>
        protected float[,] values;
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
        /// Initializes the local mask correction filter.
        /// </summary>
        public LocalCorrection()
        {
            gb = new BoxBlur(3);
            Space = Space.RGB;
            Values = new float[256, 256];
        }
        /// <summary>
        /// Initializes the local mask correction filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="values">Matrix mask</param>
        public LocalCorrection(int radius, float[,] values, Space space)
        {
            gb = new BoxBlur(radius);
            Space = space;
            Values = values;
        }
        /// <summary>
        /// Initializes the local mask correction filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="values">Matrix mask</param>
        public LocalCorrection(int width, int height, float[,] values, Space space)
        {
            gb = new BoxBlur(width, height);
            Space = space;
            Values = values;
        }
        /// <summary>
        /// Initializes the local mask correction filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="space">Color space</param>
        /// <param name="values">Matrix mask</param>
        public LocalCorrection(SizeInt size, float[,] values, Space space)
        {
            gb = new BoxBlur(size);
            Space = space;
            Values = values;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return gb.Size;
            }
            set
            {
                gb.Size = value;
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
        /// Gets or sets the matrix mask.
        /// </summary>
        public float[,] Values
        {
            get
            {
                return this.values;
            }
            set
            {
                if (value.GetLength(0) != 256 || value.GetLength(1) != 256)
                    throw new Exception("Mask size should be  256");

                this.values = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // rebuild?
            if (rebuild == true)
                this.Rebuild(); this.rebuild = false;

            // box blur
            gb.Apply(bmSrc);

            // filter
            switch (space)
            {
                case Imaging.Space.HSB:
                    ApplyHSB(bmData, bmSrc);
                    break;
                case Imaging.Space.HSL:
                    ApplyHSL(bmData, bmSrc);
                    break;
                case Imaging.Space.YCbCr:
                    ApplyYCbCr(bmData, bmSrc);
                    break;
                case Imaging.Space.RGB:
                    ApplyRGB(bmData, bmSrc);
                    break;
                default:
                    ApplyGrayscale(bmData, bmSrc);
                    break;
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
            return;
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

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyRGB(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                int i, k, k1, k2, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    p[k2] = Maths.Byte(this.values[p[k2], pSrc[k2]] * length);
                    p[k1] = Maths.Byte(this.values[p[k1], pSrc[k1]] * length);
                    p[k] = Maths.Byte(this.values[p[k], pSrc[k]] * length);
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyHSL(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                HSL lumI; HSL lumIx; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    lumI = HSL.FromRGB(p[k2], p[k1], p[k]);
                    lumIx = HSL.FromRGB(pSrc[k2], pSrc[k1], pSrc[k]);
                    lumI.Lightness = this.values[(int)(lumI.Lightness * length), (int)(lumIx.Lightness * length)];
                    rgb = lumI.ToRGB;

                    p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyHSB(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                HSB lumI; HSB lumIx; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    lumI = HSB.FromRGB(p[k2], p[k1], p[k]);
                    lumIx = HSB.FromRGB(pSrc[k2], pSrc[k1], pSrc[k]);
                    lumI.Brightness = this.values[(int)(lumI.Brightness * length), (int)(lumIx.Brightness * length)];
                    rgb = lumI.ToRGB;

                    p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyYCbCr(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                YCbCr lumI; YCbCr lumIx; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    lumI = YCbCr.FromRGB(p[k2], p[k1], p[k]);
                    lumIx = YCbCr.FromRGB(pSrc[k2], pSrc[k1], pSrc[k]);
                    lumI.Y = this.values[(int)(lumI.Y * length), (int)(lumIx.Y * length)];
                    rgb = lumI.ToRGB;

                    p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyGrayscale(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                int i, k, lumax, lumay, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    lumax = RGB.Average(p[k + 2], p[k + 1], p[k]);
                    lumay = RGB.Average(pSrc[k + 2], pSrc[k + 1], pSrc[k]);

                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(this.values[lumax, lumay] * length);
                }
            });

            return;
        }
        #endregion
    }
}
