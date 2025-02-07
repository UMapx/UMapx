using System;
using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the edge glow filter.
    /// </summary>
    [Serializable]
    public class EdgeGlow : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private readonly Erosion erosion = new Erosion();
        private readonly Operation subtraction = Operation.Subtraction;
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
            Bitmap Src0 = (Bitmap)BitmapFormat.Bitmap(bmSrc).Clone();
            BitmapData bmSrc0 = BitmapFormat.Lock32bpp(Src0);

            // Filter applying:
            erosion.Apply(bmSrc, bmSrc0);
            subtraction.Apply(bmData, bmSrc);

            // Delete resources:
            BitmapFormat.Unlock(Src0, bmSrc0);
            Src0.Dispose();
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
    }
}
