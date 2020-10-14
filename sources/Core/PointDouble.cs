using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of double numbers representing an ordered pair of X and Y coordinates.
    /// </summary>
    [Serializable]
    public struct PointDouble : ICloneable
    {
        #region Private data
        private double y;
        private double x;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of double numbers representing an ordered pair of X and Y coordinates.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        public PointDouble(double x, double y)
        {
            this.x = x;
            this.y = y;
        }
        /// <summary>
        /// Gets or sets the coordinate X.
        /// </summary>
        public double X
        {
            get
            {
                return this.x;
            }
            set
            {
                this.x = value;
            }
        }
        /// <summary>
        /// Gets or sets the coordinate Y.
        /// </summary>
        public double Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = value;
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
            return x.GetHashCode() ^ y.GetHashCode();
        }
        /// <summary>
        /// Converts a PointDouble to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", x, y);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type PointDouble.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is PointDouble) ? (this == (PointDouble)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two PointDouble objects are equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(PointDouble a, PointDouble b)
        {
            return (a.X == b.X && a.Y == b.Y);
        }
        /// <summary>
        /// Checks if two PointDouble objects are not equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(PointDouble a, PointDouble b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of PointDouble.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        object ICloneable.Clone()
        {
            return new PointDouble(x, y);
        }
        /// <summary>
        /// Creates a copy of PointDouble.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        public PointDouble Clone()
        {
            return new PointDouble(x, y);
        }
        #endregion
    }
}
