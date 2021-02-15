using System;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model YUV.
    /// </summary>
    [Serializable]
    public struct YUV : IColorSpace, ICloneable
    {
        #region Private data
        private float y;
        private float u;
        private float v;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YUV.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="u">U [-0.436, 0.436]</param>
        /// <param name="v">V [-0.614, 0.614]</param>
        public YUV(float y, float u, float v)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.u = (u > 0.436f) ? 0.436f : ((u < -0.436f) ? -0.436f : u);
            this.v = (v > 0.614f) ? 0.614f : ((v < -0.614f) ? -0.614f : v);
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
        /// Defines a component of the color model [-0.436, 0.436].
        /// </summary>
        public float U
        {
            get
            {
                return u;
            }
            set
            {
                u = (value > 0.436f) ? 0.436f : ((value < -0.436f) ? -0.436f : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.614, 0.614].
        /// </summary>
        public float V
        {
            get
            {
                return v;
            }
            set
            {
                v = (value > 0.614f) ? 0.614f : ((value < -0.614f) ? -0.614f : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YUV structure</param>
        /// <param name="item2">YUV structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YUV item1, YUV item2)
        {
            return (
                item1.Y == item2.Y
                && item1.U == item2.U
                && item1.V == item2.V
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YUV structure</param>
        /// <param name="item2">YUV structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YUV item1, YUV item2)
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

            return (this == (YUV)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ U.GetHashCode() ^ V.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + U.ToString() + "\n" + V.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YUV(this.Y, this.U, this.V);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YUV Clone()
        {
            return new YUV(this.Y, this.U, this.V);
        }
        #endregion

        #region YUV convert
        /// <summary>
        /// Converts a color model RGB in model YUV.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YUV structure</returns>
        public static YUV FromRGB(int red, int green, int blue)
        {
            float r = red / 255.0f;
            float g = green / 255.0f;
            float b = blue / 255.0f;

            float Y = 0.299f * r + 0.587f * g + 0.114f * b;
            float U = -0.147f * r - 0.288f * g + 0.436f * b;
            float V = 0.615f * r - 0.514f * g - 0.100f * b;

            return new YUV(Y, U, V);
        }
        /// <summary>
        /// Converts a color model RGB in model YUV.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YUV structure</returns>
        public static YUV FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YUV in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 1.139f * v) * 255);
                int g = (int)((y - 0.394f * u - 0.580f * v) * 255);
                int b = (int)((y + 2.032f * u) * 255);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
}
