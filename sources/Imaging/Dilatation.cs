using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the dilatation filter.
    /// </summary>
    [Serializable]
    public class Dilatation : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the dilatation filter.
        /// </summary>
        /// <param name="radius">Radius.</param>
        public Dilatation(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Initializes the dilatation filter.
        /// </summary>
        /// <param name="width">Filter width.</param>
        /// <param name="height">Filter height.</param>
        public Dilatation(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the dilatation filter.
        /// </summary>
        /// <param name="size">Filter size.</param>
        public Dilatation(SizeInt size)
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
        /// <param name="bmData">Bitmap data.</param>
        /// <param name="bmSrc">Bitmap data.</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            Morphology.Dilatation(Size.Width, Size.Height).Apply(bmData, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap.</param>
        /// <param name="Src">Bitmap.</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            try
            {
                BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
                try
                {
                    Apply(bmData, bmSrc);
                }
                finally
                {
                    BitmapFormat.Unlock(Src, bmSrc);
                }
            }
            finally
            {
                BitmapFormat.Unlock(Data, bmData);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data.</param>
        public void Apply(BitmapData bmData)
        {
            using Bitmap Src = BitmapFormat.ToBitmap(bmData);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            try
            {
                Apply(bmData, bmSrc);
            }
            finally
            {
                BitmapFormat.Unlock(Src, bmSrc);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap.</param>
        public void Apply(Bitmap Data)
        {
            using var Src = (Bitmap)Data.Clone();
            Apply(Data, Src);
        }
        #endregion
    }
}
