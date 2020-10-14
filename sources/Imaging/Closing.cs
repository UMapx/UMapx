﻿using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the closing filter.
    /// </summary>
    [Serializable]
    public class Closing : IBitmapFilter2
    {
        #region Private data
        private Erosion erosion = new Erosion();
        private Dilatation dilatation = new Dilatation();
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