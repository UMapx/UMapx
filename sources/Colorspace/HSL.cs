using System;
using UMapx.Core;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model HSL.
    /// </summary>
    [Serializable]
    public struct HSL : IColorSpace, ICloneable
    {
        #region Private data
        private float h;
        private float s;
        private float l;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure HSL.
        /// </summary>
        /// <param name="h">Hue [0, 360]</param>
        /// <param name="s">Saturation [0, 1]</param>
        /// <param name="l">Lightness [0, 1]</param>
        public HSL(float h, float s, float l)
        {
            this.h = (h > 360) ? 360 : ((h < 0) ? 0 : h);
            this.s = (s > 1) ? 1 : ((s < 0) ? 0 : s);
            this.l = (l > 1) ? 1 : ((l < 0) ? 0 : l);
        }
        /// <summary>
        /// Defines a component of the color model [0, 360].
        /// </summary>
        public float Hue
        {
            get
            {
                return h;
            }
            set
            {
                h = (value > 360) ? 360 : ((value < 0) ? 0 : value);
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
                s = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public float Lightness
        {
            get
            {
                return l;
            }
            set
            {
                l = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">HSL structure</param>
        /// <param name="item2">HSL structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(HSL item1, HSL item2)
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
        /// <param name="item1">HSL structure</param>
        /// <param name="item2">HSL structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(HSL item1, HSL item2)
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

            return (this == (HSL)obj);
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
            return new HSL(this.h, this.s, this.l);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public HSL Clone()
        {
            return new HSL(this.h, this.s, this.l);
        }
        #endregion

        #region HSL convert
        /// <summary>
        /// Converts from RGB to HSL.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>HSL structure</returns>
        public static HSL FromRGB(int red, int green, int blue)
        {
            float s = 0.0f;
            int h = 0;

            float r = red / 255.0f;
            float g = green / 255.0f;
            float b = blue / 255.0f;

            float max = Maths.Max(r, g, b);
            float min = Maths.Min(r, g, b);

            if (max == min)
            {
                h = 0;
            }
            else if (max == r && g >= b)
            {
                h = (int)(60 * (g - b) / (max - min));
            }
            else if (max == r && g < b)
            {
                h = (int)(60 * (g - b) / (max - min) + 360);
            }
            else if (max == g)
            {
                h = (int)(60 * (b - r) / (max - min) + 120);
            }
            else if (max == b)
            {
                h = (int)(60 * (r - g) / (max - min) + 240);
            }

            float l = (max + min) / 2.0f;

            if (l == 0.0 || max == min)
            {
                s = 0.0f;
            }
            else if (0.0 < l && l <= 0.5)
            {
                s = (max - min) / (max + min);
            }
            else if (l > 0.5)
            {
                s = (max - min) / (2.0f - (max + min));
            }

            return new HSL(h, s, l);
        }
        /// <summary>
        /// Converts from RGB to HSL.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>HSL structure</returns>
        public static HSL FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts from HSL to RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                if (s == 0)
                {
                    float rgb = l * 255;
                    return new RGB(rgb, rgb, rgb);
                }
                else
                {
                    float q = (l < 0.5f) ? (l * (1.0f + s)) : (l + s - (l * s));
                    float p = (2.0f * l) - q;

                    float Hk = h / 360.0f;
                    float[] T = new float[3];

                    T[0] = Hk + 0.3333f;	// Tr
                    T[1] = Hk;				// Tg
                    T[2] = Hk - 0.3333f;	// Tb

                    for (int i = 0; i < 3; i++)
                    {
                        if (T[i] < 0) T[i] += 1.0f;
                        if (T[i] > 1) T[i] -= 1.0f;

                        if ((T[i] * 6) < 1)
                        {
                            T[i] = p + ((q - p) * 6.0f * T[i]);
                        }
                        else if ((T[i] * 2.0) < 1)
                        {
                            T[i] = q;
                        }
                        else if ((T[i] * 3.0) < 2)
                        {
                            T[i] = p + (q - p) * ((2.0f / 3.0f) - T[i]) * 6.0f;
                        }
                        else T[i] = p;
                    }

                    float red = T[0] * 255;
                    float green = T[1] * 255;
                    float blue = T[2] * 255;

                    return new RGB(red, green, blue);
                }
            }
        }
        #endregion
    }
}
