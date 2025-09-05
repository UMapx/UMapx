using System;
using UMapx.Core;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model CMYK.
    /// </summary>
    [Serializable]
    public struct CMYK : IColorSpace, ICloneable
    {
        #region Private data
        private float c;
        private float m;
        private float y;
        private float k;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure CMYK.
        /// </summary>
        /// <param name="c">Cyan [0, 1]</param>
        /// <param name="m">Magenta [0, 1]</param>
        /// <param name="y">Yellow [0, 1]</param>
        /// <param name="k">Keycolor [0, 1]</param>
        public CMYK(float c, float m, float y, float k)
        {
            this.c = (c > 1) ? 1 : ((c < 0) ? 0 : c);
            this.m = (m > 1) ? 1 : ((m < 0) ? 0 : m);
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.k = (k > 1) ? 1 : ((k < 0) ? 0 : k);
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public float Cyan
        {
            get
            {
                return c;
            }
            set
            {
                c = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public float Magenta
        {
            get
            {
                return m;
            }
            set
            {
                m = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public float Yellow
        {
            get
            {
                return y;
            }
            set
            {
                y = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public float Keycolor
        {
            get
            {
                return k;
            }
            set
            {
                k = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">CMYK structure</param>
        /// <param name="item2">CMYK structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(CMYK item1, CMYK item2)
        {
            return (
                item1.Cyan == item2.Cyan
                && item1.Magenta == item2.Magenta
                && item1.Yellow == item2.Yellow
                && item1.Keycolor == item2.Keycolor
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">CMYK structure</param>
        /// <param name="item2">CMYK structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(CMYK item1, CMYK item2)
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

            return (this == (CMYK)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return c.GetHashCode() ^ m.GetHashCode() ^ y.GetHashCode() ^ k.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return $"{c}{Environment.NewLine}{m}{Environment.NewLine}{y}{Environment.NewLine}{k}";
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new CMYK(this.c, this.m, this.y, this.k);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public CMYK Clone()
        {
            return new CMYK(this.c, this.m, this.y, this.k);
        }
        #endregion

        #region CMYK convert
        /// <summary>
        /// Converts from RGB to CMYK.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>CMYK structure</returns>
        public static CMYK FromRGB(int red, int green, int blue)
        {
            float c = (255.0f - red) / 255.0f;
            float m = (255.0f - green) / 255.0f;
            float y = (255.0f - blue) / 255.0f;

            float min = Maths.Min(c, m, y);

            if (min == 1.0)
            {
                return new CMYK(0.0f, 0.0f, 0.0f, 1.0f);
            }
            else
            {
                float k = 1.0f - min;
                return new CMYK((c - min) / k, (m - min) / k, (y - min) / k, min);
            }
        }
        /// <summary>
        /// Converts from RGB to CMYK.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>CMYK structure</returns>
        public static CMYK FromRGB(RGB rgb)
        {
            return CMYK.FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts from CMYK to RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                float k1 = (1.0f - k) * 255.0f;

                int r = (int)((1.0f - c) * k1);
                int g = (int)((1.0f - m) * k1);
                int b = (int)((1.0f - y) * k1);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
}
