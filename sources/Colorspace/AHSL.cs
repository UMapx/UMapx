using System;
using UMapx.Core;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model AHSL.
    /// </summary>
    [Serializable]
    public struct AHSL : IColorSpace, ICloneable
    {
        #region Private data
        private float h;
        private float s;
        private float l;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the AHSL structure.
        /// </summary>
        /// <param name="h">Hue [0, 360]</param>
        /// <param name="s">Saturation [0, 255]</param>
        /// <param name="l">Lightness [-100, 100]</param>
        public AHSL(float h, float s, float l)
        {
            this.h = (h > 359) ? 359 : ((h < 0) ? 0 : h);
            this.s = (s > 255) ? 255 : ((s < 0) ? 0 : s);
            this.l = (l > 100) ? 100 : ((l < -100) ? -100 : l);
        }
        /// <summary>
        /// Defines a component of the color model [0, 359].
        /// </summary>
        public float Hue
        {
            get
            {
                return h;
            }
            set
            {
                h = (value > 359) ? 359 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 255].
        /// </summary>
        public float Saturation
        {
            get
            {
                return s;
            }
            set
            {
                s = (value > 255) ? 255 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-100, 100].
        /// </summary>
        public float Lightness
        {
            get
            {
                return l;
            }
            set
            {
                l = (value > 100) ? 100 : ((value < -100) ? -100 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">AHSL structure</param>
        /// <param name="item2">AHSL structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(AHSL item1, AHSL item2)
        {
            return (
                item1.Hue == item2.Hue
                && item1.Saturation == item2.Saturation
                && item1.Lightness == item2.Lightness
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">AHSL structure</param>
        /// <param name="item2">AHSL structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(AHSL item1, AHSL item2)
        {
            return !(item1 == item2);
        }
        #endregion

        #region Methods
        /// <summary>
        /// Defines whether the specified System.Object is equal to the current System.Object.
        /// </summary>
        /// <param name="obj">Element</param>
        /// <returns>Boolean</returns>
        public override bool Equals(Object obj)
        {
            if (obj == null || GetType() != obj.GetType()) return false;

            return (this == (AHSL)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return h.GetHashCode() ^ s.GetHashCode() ^ l.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return $"{h}{Environment.NewLine}{s}{Environment.NewLine}{l}";
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new AHSL(this.h, this.s, this.l);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public AHSL Clone()
        {
            return new AHSL(this.h, this.s, this.l);
        }
        #endregion

        #region AHSL convert
        /// <summary>
        /// Converts a color model RGB in model AHSL.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>AHSL structure</returns>
        public static AHSL FromRGB(int red, int green, int blue)
        {
            float s = 0, l = 0, h = 0;

            float max = Maths.Max(red, green, blue);
            float min = Maths.Min(red, green, blue);

            if (max == min)
            {
                h = 0;
            }
            else if (max == red && green >= blue)
            {
                h = 60 * (green - blue) / (max - min);
            }
            else if (max == red && green < blue)
            {
                h = 60 * (green - blue) / (max - min) + 360;
            }
            else if (max == green)
            {
                h = 60 * (blue - red) / (max - min) + 120;
            }
            else if (max == blue)
            {
                h = 60 * (red - green) / (max - min) + 240;
            }

            float gray = (red + green + blue) / 3.0f;
            float r0 = 0, g0 = 0, b0 = 0;

            if ((h >= 0) && (h <= 60)) { r0 = 255.0f; g0 = 4.25f * h; }
            else if ((h > 60) && (h <= 120)) { g0 = 255.0f; r0 = 255 - 4.25f * (h - 60); }
            else if ((h > 120) && (h <= 180)) { g0 = 255.0f; b0 = 4.25f * (h - 120); }
            else if ((h > 180) && (h <= 240)) { b0 = 255.0f; g0 = 255 - 4.25f * (h - 180); }
            else if ((h > 240) && (h <= 300)) { b0 = 255.0f; r0 = 4.25f * (h - 240); }
            else if ((h > 300) && (h <= 360)) { r0 = 255.0f; b0 = 255 - 4.25f * (h - 300); }

            float gray0 = (r0 + g0 + b0) / 3.0f;

            if (gray == gray0) { l = 0; }
            else if (gray > gray0) { l = 100 * (gray - gray0) / (255.0f - gray0); }
            else if (gray < gray0) { l = 100 * (gray - gray0) / gray0; }

            if (l > 0) { r0 = r0 + l * (255 - r0) / 100.0f; }
            else if (l < 0) { r0 = r0 + l * r0 / 100.0f; }

            if (red == gray) { s = 0; }
            else { s = 255 * Math.Abs(red - gray) / (Math.Abs(r0 - gray)); }

            return new AHSL(h, s, l);
        }
        /// <summary>
        /// Converts a color model RGB in model AHSL.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>AHSL structure</returns>
        public static AHSL FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model AHSL in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                float r = 0, g = 0, b = 0;

                if ((h >= 0) && (h <= 60)) { r = 255.0f; g = 4.25f * h; }
                else if ((h > 60) && (h <= 120)) { g = 255.0f; r = 255 - 4.25f * (h - 60); }
                else if ((h > 120) && (h <= 180)) { g = 255.0f; b = 4.25f * (h - 120); }
                else if ((h > 180) && (h <= 240)) { b = 255.0f; g = 255 - 4.25f * (h - 180); }
                else if ((h > 240) && (h <= 300)) { b = 255.0f; r = 4.25f * (h - 240); }
                else if ((h > 300) && (h <= 360)) { r = 255.0f; b = 255 - 4.25f * (h - 300); }

                if (l > 0)
                {
                    r = r + l * (255 - r) / 100;
                    g = g + l * (255 - g) / 100;
                    b = b + l * (255 - b) / 100;
                }
                else
                {
                    r = r + l * r / 100.0f;
                    g = g + l * g / 100.0f;
                    b = b + l * b / 100.0f;
                }

                float gray = (r + g + b) / 3.0f;
                r = gray + s * (r - gray) / 255.0f;
                g = gray + s * (g - gray) / 255.0f;
                b = gray + s * (b - gray) / 255.0f;

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
}
