using System;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model YPbPr.
    /// </summary>
    [Serializable]
    public struct YPbPr : IColorSpace, ICloneable
    {
        #region Private data
        private float y;
        private float pb;
        private float pr;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YPbPr.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="pb">Pb [-0.5, 0.5]</param>
        /// <param name="pr">Pr [-0.5, 0.5]</param>
        public YPbPr(float y, float pb, float pr)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.pb = (pb > 0.5) ? 0.5f : ((pb < -0.5) ? -0.5f : pb);
            this.pr = (pr > 0.5) ? 0.5f : ((pr < -0.5) ? -0.5f : pr);
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
        public float Pb
        {
            get
            {
                return pb;
            }
            set
            {
                pb = (value > 0.5) ? 0.5f : ((value < -0.5) ? -0.5f : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.5, 0.5].
        /// </summary>
        public float Pr
        {
            get
            {
                return pr;
            }
            set
            {
                pr = (value > 0.5) ? 0.5f : ((value < -0.5) ? -0.5f : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YPbPr structure</param>
        /// <param name="item2">YPbPr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YPbPr item1, YPbPr item2)
        {
            return (
                item1.Y == item2.Y
                && item1.Pb == item2.Pb
                && item1.Pr == item2.Pr
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YPbPr structure</param>
        /// <param name="item2">YPbPr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YPbPr item1, YPbPr item2)
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

            return (this == (YPbPr)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ Pb.GetHashCode() ^ Pr.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + Pb.ToString() + "\n" + Pr.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YPbPr(this.Y, this.Pb, this.Pr);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YPbPr Clone()
        {
            return new YPbPr(this.Y, this.Pb, this.Pr);
        }
        #endregion

        #region YPbPr convert
        /// <summary>
        /// Converts a color model RGB in model YPbPr.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YPbPr structure</returns>
        public static YPbPr FromRGB(int red, int green, int blue)
        {
            float r = red / 255.0f;
            float g = green / 255.0f;
            float b = blue / 255.0f;

            float Y = 0.299f * r + 0.587f * g + 0.114f * b;
            float Cb = -0.169f * r - 0.331f * g + 0.500f * b;
            float Cr = 0.500f * r - 0.419f * g - 0.081f * b;

            return new YPbPr(Y, Cb, Cr);
        }
        /// <summary>
        /// Converts a color model RGB in model YPbPr.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YPbPr structure</returns>
        public static YPbPr FromRGB(RGB rgb)
        {
            return YPbPr.FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YPbPr in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 1.402 * pr) * 255.0);
                int g = (int)((y - 0.344 * pb - 0.714 * pr) * 255.0);
                int b = (int)((y + 1.772 * pb) * 255.0);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
}
