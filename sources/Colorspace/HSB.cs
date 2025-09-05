using System;
using UMapx.Core;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model HSB.
    /// </summary>
    [Serializable]
    public struct HSB : IColorSpace, ICloneable
    {
        #region Private data
        private float h;
        private float s;
        private float b;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure HSB.
        /// </summary>
        /// <param name="h">Hue [0, 359]</param>
        /// <param name="s">Saturation [0, 1]</param>
        /// <param name="b">Brightness [0, 1]</param>
        public HSB(float h, float s, float b)
        {
            this.h = (h > 359) ? 359 : ((h < 0) ? 0 : h);
            this.s = (s > 1) ? 1 : ((s < 0) ? 0 : s);
            this.b = (b > 1) ? 1 : ((b < 0) ? 0 : b);
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
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public float Saturation
        {
            get
            {
                return s;
            }
            set
            {
                s = (value > 1) ? 1 : ((value < 0.00001f) ? 0.00001f : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public float Brightness
        {
            get
            {
                return b;
            }
            set
            {
                b = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">HSB structure</param>
        /// <param name="item2">HSB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(HSB item1, HSB item2)
        {
            return (
                item1.Hue == item2.Hue
                && item1.Saturation == item2.Saturation
                && item1.Brightness == item2.Brightness
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">HSB structure</param>
        /// <param name="item2">HSB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(HSB item1, HSB item2)
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

            return (this == (HSB)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return h.GetHashCode() ^ s.GetHashCode() ^ b.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return $"{h}{Environment.NewLine}{s}{Environment.NewLine}{b}";
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new HSB(this.h, this.s, this.b);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public HSB Clone()
        {
            return new HSB(this.h, this.s, this.b);
        }
        #endregion

        #region HSB convert
        /// <summary>
        /// Converts from RGB to HSB.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>HSB structure</returns>
        public static HSB FromRGB(int red, int green, int blue)
        {
            float r = red / 255.0f;
            float g = green / 255.0f;
            float b = blue / 255.0f;

            float max = Maths.Max(r, g, b);
            float min = Maths.Min(r, g, b);

            float h = 0;
            float l = max - min;

            if (max == r && g >= b)
            {
                if (l == 0.0) h = 0;
                else h = (int)(60 * (g - b) / l);
            }
            else if (max == r && g < b)
            {
                h = (int)(60 * (g - b) / l + 360);
            }
            else if (max == g)
            {
                h = (int)(60 * (b - r) / l + 120);
            }
            else if (max == b)
            {
                h = (int)(60 * (r - g) / l + 240);
            }

            float s = (max == 0.0f) ? 0.0f : (1 - (min / max));
            return new HSB(h, s, max);
        }
        /// <summary>
        /// Converts from RGB to HSB.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>HSB structure</returns>
        public static HSB FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts from HSB to RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                float red = 0;
                float green = 0;
                float blue = 0;

                if (b <= 0)
                {
                    red = green = blue = 0;
                }

                if (s <= 0)
                {
                    red = green = blue = b;
                }
                else
                {
                    float Hi = h / 60.0f;
                    int HiMod6 = (int)(Math.Floor(Hi));
                    float Sector = Hi - HiMod6;

                    float Vmin = b * (1 - s);
                    float Vdec = b * (1 - (s * Sector));
                    float Vinc = b * (1 - (s * (1 - Sector)));

                    switch (HiMod6)
                    {
                        case 0:
                            red = b;
                            green = Vinc;
                            blue = Vmin;
                            break;
                        case 1:
                            red = Vdec;
                            green = b;
                            blue = Vmin;
                            break;
                        case 2:
                            red = Vmin;
                            green = b;
                            blue = Vinc;
                            break;
                        case 3:
                            red = Vmin;
                            green = Vdec;
                            blue = b;
                            break;
                        case 4:
                            red = Vinc;
                            green = Vmin;
                            blue = b;
                            break;
                        case 5:
                            red = b;
                            green = Vmin;
                            blue = Vdec;
                            break;
                    }
                }
                return new RGB(red * 255, green * 255, blue * 255);
            }
        }
        #endregion
    }
}
