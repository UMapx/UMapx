using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of double numbers representing an ordered pair of width and height.
    /// </summary>
    [Serializable]
    public struct SizeDouble : ICloneable
    {
        #region Private data
        private double width;
        private double height;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of double numbers representing an ordered pair of width and height.
        /// </summary>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        public SizeDouble(double width, double height)
        {
            this.height = height;
            this.width = width;
        }
        /// <summary>
        /// Gets or sets the height.
        /// </summary>
        public double Height
        {
            get
            {
                return this.height;
            }
            set
            {
                this.height = value;
            }
        }
        /// <summary>
        /// Gets or sets the width.
        /// </summary>
        public double Width
        {
            get
            {
                return this.width;
            }
            set
            {
                this.width = value;
            }
        }
        #endregion

        #region Overrides
        /// <summary>
        /// Returns the hash code for this object.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return width.GetHashCode() ^ height.GetHashCode();
        }
        /// <summary>
        /// Converts a SizeDouble to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", width, height);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type SizeDouble.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is SizeDouble) ? (this == (SizeDouble)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two SizeDouble objects are equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(SizeDouble a, SizeDouble b)
        {
            return (a.Width == b.Width && a.Height == b.Height);
        }
        /// <summary>
        /// Checks if two SizeDouble objects are not equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(SizeDouble a, SizeDouble b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of SizeDouble.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        object ICloneable.Clone()
        {
            return new SizeDouble(width, height);
        }
        /// <summary>
        /// Creates a copy of SizeDouble.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        public SizeDouble Clone()
        {
            return new SizeDouble(width, height);
        }
        #endregion
    }
}
