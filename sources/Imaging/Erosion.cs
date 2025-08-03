using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the erosion filter.
    /// </summary>
    [Serializable]
    public class Erosion : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public Erosion(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public Erosion(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public Erosion(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(rw, rh);
            }
            set
            {
                this.rw = value.Width;
                this.rh = value.Height;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            Morphology.Erosion(Size.Width, Size.Height).Apply(bmData, bmSrc);
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
        }
        #endregion
    }
}
