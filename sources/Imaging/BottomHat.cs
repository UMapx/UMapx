using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the bottom-hat filter.
    /// </summary>
    [Serializable]
    public class BottomHat : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private readonly Closing closing = new Closing();
        private readonly Operation subtraction = Operation.Subtraction;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the bottom-hat filter.
        /// </summary>
        /// <param name="radius">Radius.</param>
        public BottomHat(int radius = 3)
        {
            closing = new Closing(radius);
        }
        /// <summary>
        /// Initializes the bottom-hat filter.
        /// </summary>
        /// <param name="width">Filter width.</param>
        /// <param name="height">Filter height.</param>
        public BottomHat(int width, int height)
        {
            closing = new Closing(width, height);
        }
        /// <summary>
        /// Initializes the bottom-hat filter.
        /// </summary>
        /// <param name="size">Filter size.</param>
        public BottomHat(SizeInt size)
        {
            closing = new Closing(size);
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return closing.Size;
            }
            set
            {
                closing.Size = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data.</param>
        /// <param name="bmSrc">Bitmap data.</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // Creating resources:
            using Bitmap Src0 = BitmapFormat.ToBitmap(bmSrc);
            BitmapData bmSrc0 = BitmapFormat.Lock32bpp(Src0);
            try
            {
                // Filter applying:
                closing.Apply(bmSrc, bmSrc0);
                subtraction.Apply(bmData, bmSrc);
            }
            finally
            {
                BitmapFormat.Unlock(Src0, bmSrc0);
            }
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
