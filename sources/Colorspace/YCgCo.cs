using System;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model YCgCo.
    /// </summary>
    [Serializable]
    public struct YCgCo : IColorSpace, ICloneable
    {
        #region Private data
        private float y;
        private float cg;
        private float co;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YCgCo.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="cg">Cg [-0.5, 0.5]</param>
        /// <param name="co">Co [-0.5, 0.5]</param>
        public YCgCo(float y, float cg, float co)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.cg = (cg > 0.5) ? 0.5f : ((cg < -0.5) ? -0.5f : cg);
            this.co = (co > 0.5) ? 0.5f : ((co < -0.5) ? -0.5f : co);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public float Y
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
        /// Defines a component of the color model [-0.5, 0.5].
        /// </summary>
        public float Cg
        {
            get
            {
                return cg;
            }
            set
            {
                cg = (value > 0.5) ? 0.5f : ((value < -0.5) ? -0.5f : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.5, 0.5].
        /// </summary>
        public float Co
        {
            get
            {
                return co;
            }
            set
            {
                co = (value > 0.5) ? 0.5f : ((value < -0.5) ? -0.5f : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YCgCo structure</param>
        /// <param name="item2">YCgCo structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YCgCo item1, YCgCo item2)
        {
            return (
                item1.Y == item2.Y
                && item1.Cg == item2.Cg
                && item1.Co == item2.Co
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YCgCo structure</param>
        /// <param name="item2">YCgCo structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YCgCo item1, YCgCo item2)
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

            return (this == (YCgCo)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return y.GetHashCode() ^ cg.GetHashCode() ^ co.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return $"{y}{Environment.NewLine}{cg}{Environment.NewLine}{co}";
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YCgCo(this.y, this.cg, this.co);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YCgCo Clone()
        {
            return new YCgCo(this.y, this.cg, this.co);
        }
        #endregion

        #region YCgCo convert
        /// <summary>
        /// Converts from RGB to YCgCo.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YCgCo structure</returns>
        public static YCgCo FromRGB(int red, int green, int blue)
        {
            float r = red / 255.0f;
            float g = green / 255.0f;
            float b = blue / 255.0f;

            float Y = 0.25f * r + 0.5f * g + 0.25f * b;
            float Cg = -0.25f * r + 0.5f * g - 0.25f * b;
            float Co = 0.5f * r - 0.0f * g - 0.5f * b;

            return new YCgCo(Y, Cg, Co);
        }
        /// <summary>
        /// Converts from RGB to YCgCo.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YCgCo structure</returns>
        public static YCgCo FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts from YCgCo to RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y - cg + co) * 255.0);
                int g = (int)((y + cg) * 255.0);
                int b = (int)((y - cg - co) * 255.0);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
}
