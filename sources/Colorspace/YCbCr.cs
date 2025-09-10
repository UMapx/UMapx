using System;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model YCbCr.
    /// </summary>
    [Serializable]
    public struct YCbCr : IColorSpace, ICloneable
    {
        #region Private data
        private float y;
        private float cb;
        private float cr;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YCbCr.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="cb">Cb [0, 1]</param>
        /// <param name="cr">Cr [0, 1]</param>
        public YCbCr(float y, float cb, float cr)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.cb = (cb > 1) ? 1 : ((cb < 0) ? 0 : cb);
            this.cr = (cr > 1) ? 1 : ((cr < 0) ? 0 : cr);
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
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public float Cb
        {
            get
            {
                return cb;
            }
            set
            {
                cb = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public float Cr
        {
            get
            {
                return cr;
            }
            set
            {
                cr = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YCbCr structure</param>
        /// <param name="item2">YCbCr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YCbCr item1, YCbCr item2)
        {
            return (
                item1.Y == item2.Y
                && item1.Cb == item2.Cb
                && item1.Cr == item2.Cr
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YCbCr structure</param>
        /// <param name="item2">YCbCr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YCbCr item1, YCbCr item2)
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

            return (this == (YCbCr)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return y.GetHashCode() ^ cb.GetHashCode() ^ cr.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return $"{y}{Environment.NewLine}{cb}{Environment.NewLine}{cr}";
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YCbCr(this.y, this.cb, this.cr);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YCbCr Clone()
        {
            return new YCbCr(this.y, this.cb, this.cr);
        }
        #endregion

        #region YCbCr convert
        /// <summary>
        /// Converts from RGB to YCbCr.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YCbCr structure</returns>
        public static YCbCr FromRGB(int red, int green, int blue)
        {
            float r = red / 255.0f;
            float g = green / 255.0f;
            float b = blue / 255.0f;

            float Y = 0.299f * r + 0.587f * g + 0.114f * b;
            float Cb = -0.168736f * r - 0.331264f * g + 0.5f * b + 0.5f;
            float Cr = 0.5f * r - 0.418688f * g - 0.081312f * b + 0.5f;

            return new YCbCr(Y, Cb, Cr);
        }
        /// <summary>
        /// Converts from RGB to YCbCr.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YCbCr structure</returns>
        public static YCbCr FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts from YCbCr to RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 1.402f * (cr - 0.5f)) * 255);
                int g = (int)((y - 0.344136f * (cb - 0.5f) - 0.714136f * (cr - 0.5f)) * 255);
                int b = (int)((y + 1.772f * (cb - 0.5f)) * 255);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
}
