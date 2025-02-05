using System;
using UMapx.Core;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the closing filter.
    /// </summary>
    [Serializable]
    public class Closing : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private readonly Erosion erosion = new Erosion();
        private readonly Dilatation dilatation = new Dilatation();
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the closing filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public Closing(int radius = 3)
        {
            erosion = new Erosion(radius);
            dilatation = new Dilatation(radius);
        }
        /// <summary>
        /// Initializes the closing filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public Closing(int width, int height)
        {
            erosion = new Erosion(width, height);
            dilatation = new Dilatation(width, height);
        }
        /// <summary>
        /// Initializes the closing filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public Closing(SizeInt size)
        {
            erosion = new Erosion(size);
            dilatation = new Dilatation(size);
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return erosion.Size;
            }
            set
            {
                erosion.Size = value;
                dilatation.Size = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            dilatation.Apply(bmSrc, bmData);
            erosion.Apply(bmData, bmSrc);
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
    }
}
