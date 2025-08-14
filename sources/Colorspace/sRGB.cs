using System;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model sRGB.
    /// </summary>
    [Serializable]
    public struct sRGB : IColorSpace, ICloneable
    {
        #region private data
        private float r;
        private float g;
        private float b;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure sRGB.
        /// </summary>
        /// <param name="red">Red [0, 1]</param>
        /// <param name="green">Green [0, 1]</param>
        /// <param name="blue">Blue [0, 1]</param>
        public sRGB(float red, float green, float blue)
        {
            this.r = (red > 1) ? 1 : ((red < 0) ? 0 : red);
            this.g = (green > 1) ? 1 : ((green < 0) ? 0 : green);
            this.b = (blue > 1) ? 1 : ((blue < 0) ? 0 : blue);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public float Red
        {
            get
            {
                return r;
            }
            set
            {
                r = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public float Green
        {
            get
            {
                return g;
            }
            set
            {
                g = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public float Blue
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
        /// <param name="item1">sRGB structure</param>
        /// <param name="item2">sRGB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(sRGB item1, sRGB item2)
        {
            return (
                item1.Red == item2.Red
                && item1.Green == item2.Green
                && item1.Blue == item2.Blue
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">sRGB structure</param>
        /// <param name="item2">sRGB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(sRGB item1, sRGB item2)
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

            return (this == (sRGB)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return r.GetHashCode() ^ g.GetHashCode() ^ b.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return $"{r}{Environment.NewLine}{g}{Environment.NewLine}{b}";
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new sRGB(this.r, this.g, this.b);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public sRGB Clone()
        {
            return new sRGB(this.r, this.g, this.b);
        }
        #endregion

        #region sRGB convert
        /// <summary>
        /// Converts a color model RGB in model sRGB.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>sRGB structure</returns>
        public static sRGB FromRGB(int red, int green, int blue)
        {
            float r = red / 255.0f;
            float g = green / 255.0f;
            float b = blue / 255.0f;

            return new sRGB(r, g, b);
        }
        /// <summary>
        /// Converts a color model RGB in model sRGB.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>sRGB structure</returns>
        public static sRGB FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model sRGB in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                return new RGB(r * 255.0f, g * 255.0f, b * 255.0f);

            }
        }
        #endregion
    }
}
