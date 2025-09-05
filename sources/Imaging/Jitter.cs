﻿using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the jitter filter.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    [Serializable]
    public class Jitter : PointAddition, IBitmapFilter2
    {
        #region Private data
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the jitter filter.
        /// </summary>
        /// <param name="value">Value [0, 100]</param>
        public Jitter(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Initializes the jitter filter.
        /// </summary>
        public Jitter() { }
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
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Noise(this.width, this.height, this.value);
        }
        #endregion
    }
}
