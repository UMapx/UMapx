﻿using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the brightness correction filter.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://esate.ru/uroki/OpenGL/image_processing/_p4106/
    /// </remarks>
    [Serializable]
    public class BrightnessCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private float brightness;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the brightness correction filter.
        /// </summary>
        /// <param name="brightness">Brightness [-1, 1]</param>
        /// <param name="space">Color space</param>
        public BrightnessCorrection(float brightness, Space space)
        {
            Brightness = brightness; this.Space = space;
        }
        /// <summary>
        /// Initializes the brightness correction filter.
        /// </summary>
        public BrightnessCorrection()
        {
            Brightness = 0.5f;
        }
        /// <summary>
        /// Gets or sets the brightness value [-1, 1].
        /// </summary>
        public float Brightness
        {
            get
            {
                return this.brightness;
            }
            set
            {
                this.brightness = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Add(this.brightness / 2.0f, 256);
        }
        #endregion
    }
}
