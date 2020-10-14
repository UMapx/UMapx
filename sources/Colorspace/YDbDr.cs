using System;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model YDbDr.
    /// </summary>
    [Serializable]
    public struct YDbDr : IColorSpace, ICloneable
    {
        #region Private data
        private double y;
        private double db;
        private double dr;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YDbDr.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="db">Db [-1.333, 1.333]</param>
        /// <param name="dr">Dr [-1.333, 1.333]</param>
        public YDbDr(double y, double db, double dr)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.db = (db > 1.333) ? 1.333 : ((db < -1.333) ? -1.333 : db);
            this.dr = (dr > 1.333) ? 1.333 : ((dr < -1.333) ? -1.333 : dr);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Y
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
        /// Defines a component of the color model [-1.333, 1.333].
        /// </summary>
        public double Db
        {
            get
            {
                return db;
            }
            set
            {
                db = (value > 1.333) ? 1 : ((value < -1.333) ? -1.333 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-1.333, 1.333].
        /// </summary>
        public double Dr
        {
            get
            {
                return dr;
            }
            set
            {
                dr = (value > 1.333) ? 1.333 : ((value < -1.333) ? -1.333 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YDbDr structure</param>
        /// <param name="item2">YDbDr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YDbDr item1, YDbDr item2)
        {
            return (
                item1.Y == item2.Y
                && item1.Db == item2.Db
                && item1.Dr == item2.Dr
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YDbDr structure</param>
        /// <param name="item2">YDbDr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YDbDr item1, YDbDr item2)
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

            return (this == (YDbDr)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ Db.GetHashCode() ^ Dr.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + Db.ToString() + "\n" + Dr.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YDbDr(this.Y, this.Db, this.Dr);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YDbDr Clone()
        {
            return new YDbDr(this.Y, this.Db, this.Dr);
        }
        #endregion

        #region YDbDr convert
        /// <summary>
        /// Converts a color model RGB in model YDbDr.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YDbDr structure</returns>
        public static YDbDr FromRGB(int red, int green, int blue)
        {
            double r = red / 255.0f;
            double g = green / 255.0f;
            double b = blue / 255.0f;

            double Y = 0.299f * r + 0.587f * g + 0.114f * b;
            double Cb = -0.450f * r - 0.883f * g + 1.333f * b;
            double Cr = -1.333f * r - 1.116f * g - 0.217f * b;

            return new YDbDr(Y, Cb, Cr);
        }
        /// <summary>
        /// Converts a color model RGB in model YDbDr.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YDbDr structure</returns>
        public static YDbDr FromRGB(RGB rgb)
        {
            return YDbDr.FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YDbDr in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 0.00009 * db - 0.52591 * dr) * 255);
                int g = (int)((y - 0.12913 * db + 0.26789 * dr) * 255);
                int b = (int)((y + 0.66467 * db - 0.00007 * dr) * 255);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
}
