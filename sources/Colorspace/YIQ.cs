using System;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model YIQ.
    /// </summary>
    [Serializable]
    public struct YIQ : IColorSpace, ICloneable
    {
        #region Private data
        private float y;
        private float i;
        private float q;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YIQ.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="i">I [-0.5957, 0.5957]</param>
        /// <param name="q">Q [-0.5226, 0.5226]</param>
        public YIQ(float y, float i, float q)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.i = (i > 0.5957f) ? 0.5957f : ((i < -0.5957f) ? -0.5957f : i);
            this.q = (q > 0.5226f) ? 0.5226f : ((q < -0.5226f) ? -0.5226f : q);
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
        /// Defines a component of the color model [-0.5957, 0.5957].
        /// </summary>
        public float I
        {
            get
            {
                return i;
            }
            set
            {
                i = (value > 0.5957f) ? 0.5957f : ((value < -0.5957f) ? -0.5957f : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.5226, 0.5226].
        /// </summary>
        public float Q
        {
            get
            {
                return q;
            }
            set
            {
                q = (value > 0.5226f) ? 0.5226f : ((value < -0.5226f) ? -0.5226f : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YIQ structure</param>
        /// <param name="item2">YIQ structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YIQ item1, YIQ item2)
        {
            return (
                item1.Y == item2.Y
                && item1.I == item2.I
                && item1.Q == item2.Q
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YIQ structure</param>
        /// <param name="item2">YIQ structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YIQ item1, YIQ item2)
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

            return (this == (YIQ)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ I.GetHashCode() ^ Q.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + I.ToString() + "\n" + Q.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YIQ(this.Y, this.I, this.Q);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YIQ Clone()
        {
            return new YIQ(this.Y, this.I, this.Q);
        }
        #endregion

        #region YIQ convert
        /// <summary>
        /// Converts a color model RGB in model YIQ.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YIQ structure</returns>
        public static YIQ FromRGB(int red, int green, int blue)
        {
            float r = red / 255.0f;
            float g = green / 255.0f;
            float b = blue / 255.0f;

            float Y = 0.299f * r + 0.587f * g + 0.114f * b;
            float I = 0.596f * r - 0.274f * g - 0.322f * b;
            float Q = 0.211f * r - 0.522f * g + 0.311f * b;

            return new YIQ(Y, I, Q);
        }
        /// <summary>
        /// Converts a color model RGB in model YIQ.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YIQ structure</returns>
        public static YIQ FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YIQ in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 0.956f * i + 0.623f * q) * 255);
                int g = (int)((y - 0.272f * i - 0.648f * q) * 255);
                int b = (int)((y - 1.105f * i + 1.705f * q) * 255);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
}
