using System;
using UMapx.Core;
using SkiaDrawing;

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
        /// <param name="radius">Radius</param>
        public BottomHat(int radius = 3)
        {
            closing = new Closing(radius);
        }
        /// <summary>
        /// Initializes the bottom-hat filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public BottomHat(int width, int height)
        {
            closing = new Closing(width, height);
        }
        /// <summary>
        /// Initializes the bottom-hat filter.
        /// </summary>
        /// <param name="size">Filter size</param>
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
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // Creating resources:
            Bitmap Src0 = (Bitmap)BitmapFormat.Bitmap(bmSrc).Clone();
            BitmapData bmSrc0 = BitmapFormat.Lock32bpp(Src0);

            // Filter applying:
            closing.Apply(bmSrc, bmSrc0);
            subtraction.Apply(bmData, bmSrc);

            // Delete resources:
            BitmapFormat.Unlock(Src0, bmSrc0);
            Src0.Dispose();
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
