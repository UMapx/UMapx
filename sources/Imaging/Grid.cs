﻿using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the grid filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Grid : PointAddition, IBitmapFilter2
    {
        #region Private data
        private int value = 20;
        private int thickness = 1;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the grid filter.
        /// </summary>
        /// <param name="value">Value [0, 100]</param>
        /// <param name="thickness">Thickness (>0)</param>
        public Grid(int value, int thickness = 1)
        {
            Value = value;
            Thickness = thickness;
        }
        /// <summary>
        /// Initializes the grid filter.
        /// </summary>
        public Grid() { }
        /// <summary>
        /// Gets or sets the value.
        /// </summary>
        public int Value
        {
            get
            {
                return this.value;
            }
            set
            {
                this.value = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the thickness (>0).
        /// </summary>
        public int Thickness
        {
            get
            {
                return this.thickness;
            }
            set
            {
                this.thickness = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Grid(this.width, this.height, this.value, this.thickness);
        }
        #endregion
    }
}
