using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the edge glow filter.
    /// </summary>
    [Serializable]
    public class EdgeGlow : IBitmapFilter2
    {
        #region Private data
        private Erosion erosion = new Erosion();
        private Operation subtraction = Operation.Subtraction;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the edge glow filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public EdgeGlow(int radius = 3)
        {
            erosion = new Erosion(radius);
        }
        /// <summary>
        /// Initializes the edge glow filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public EdgeGlow(int width, int height)
        {
            erosion = new Erosion(width, height);
        }
        /// <summary>
        /// Initializes the edge glow filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public EdgeGlow(SizeInt size)
        {
            erosion = new Erosion(size);
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return this.erosion.Size;
            }
            set
            {
                this.erosion.Size = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // Creating resources:
            Bitmap Src0 = (Bitmap)BitmapConverter.Bitmap(bmSrc).Clone();
            BitmapData bmSrc0 = BitmapConverter.Lock32bpp(Src0);

            // Filter applying:
            erosion.Apply(bmSrc, bmSrc0);
            subtraction.Apply(bmData, bmSrc);

            // Delete resources:
            BitmapConverter.Unlock(Src0, bmSrc0);
            Src0.Dispose();
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
}
