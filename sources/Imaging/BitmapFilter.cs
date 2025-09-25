using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Transform;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the bitmap filter.
    /// </summary>
    [Serializable]
    public class BitmapFilter : IBitmapFilter
    {
        #region Private data
        private IFilter filter;
        private Space space;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the bitmap filter.
        /// </summary>
        /// <param name="filter">Filter</param>
        /// <param name="space">Color space</param>
        public BitmapFilter(IFilter filter, Space space = Space.RGB)
        {
            this.filter = filter;
            this.space = space;
        }
        /// <summary>
        /// Gets or sets the filter.
        /// </summary>
        public IFilter Filter
        {
            get
            {
                return this.filter;
            }
            set
            {
                this.filter = value;
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
        public void Apply(BitmapData bmData)
        {
            if (bmData.PixelFormat != PixelFormat.Format32bppArgb)
                throw new NotSupportedException("Only support Format32bppArgb pixelFormat");

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

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private void ApplyRGB(BitmapData bmData)
        {
            float[][,] rgb = BitmapMatrix.ToRGB(bmData, true);

            this.filter.Apply(rgb[0]);
            this.filter.Apply(rgb[1]);
            this.filter.Apply(rgb[2]);

            BitmapMatrix.FromRGB(rgb, bmData);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private void ApplyHSB(BitmapData bmData)
        {
            float[][,] hsb = BitmapMatrix.ToHSB(bmData, true);
            this.filter.Apply(hsb[2]);
            BitmapMatrix.FromHSB(hsb, bmData);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private void ApplyHSL(BitmapData bmData)
        {
            float[][,] hsl = BitmapMatrix.ToHSL(bmData, true);
            this.filter.Apply(hsl[2]);
            BitmapMatrix.FromHSL(hsl, bmData);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private void ApplyYCbCr(BitmapData bmData)
        {
            float[][,] ycbcr = BitmapMatrix.ToYCbCr(bmData, true);
            this.filter.Apply(ycbcr[0]);
            BitmapMatrix.FromYCbCr(ycbcr, bmData);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private void ApplyGrayscale(BitmapData bmData)
        {
            float[,] y = BitmapMatrix.ToGrayscale(bmData);
            this.filter.Apply(y);
            BitmapMatrix.FromGrayscale(y, bmData);
        }
        #endregion
    }
}
