using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of float numbers representing an ordered pair of width and height.
    /// </summary>
    [Serializable]
    public struct SizeFloat : ICloneable
    {
        #region Private data
        private float width;
        private float height;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of float numbers representing an ordered pair of width and height.
        /// </summary>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        public SizeFloat(float width, float height)
        {
            this.height = height;
            this.width = width;
        }
        /// <summary>
        /// Gets or sets the height.
        /// </summary>
        public float Height
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
        public float Width
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
        /// Converts a SizeFloat to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", width, height);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type SizeFloat.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is SizeFloat) ? (this == (SizeFloat)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two SizeFloat objects are equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(SizeFloat a, SizeFloat b)
        {
            return (a.Width == b.Width && a.Height == b.Height);
        }
        /// <summary>
        /// Checks if two SizeFloat objects are not equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(SizeFloat a, SizeFloat b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of SizeFloat.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        object ICloneable.Clone()
        {
            return new SizeFloat(width, height);
        }
        /// <summary>
        /// Creates a copy of SizeFloat.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        public SizeFloat Clone()
        {
            return new SizeFloat(width, height);
        }
        #endregion
    }
}
