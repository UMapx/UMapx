using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of double numbers representing a line segment.
    /// </summary>
    [Serializable]
    public struct RangeDouble : ICloneable
    {
        #region Private data
        private double max;
        private double min;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of double numbers representing a line segment.
        /// </summary>
        /// <param name="min">Lower bound of the segment</param>
        /// <param name="max">Upper bound of the segment</param>
        public RangeDouble(double min, double max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Gets or sets the lower bound of the line segment.
        /// </summary>
        public double Min
        {
            get
            {
                return this.min;
            }
            set
            {
                this.min = value;
            }
        }
        /// <summary>
        /// Gets or sets the upper bound of the line segment.
        /// </summary>
        public double Max
        {
            get
            {
                return this.max;
            }
            set
            {
                this.max = value;
            }
        }
        /// <summary>
        /// Checks if the value is in the specified interval.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Boolean</returns>
        public bool IsOnRange(double x)
        {
            if ((x >= this.min) && (x <= this.max))
            {
                return true;
            }
            return false;
        }
        #endregion

        #region Overrides
        /// <summary>
        /// Returns the hash code for this object.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return min.GetHashCode() ^ max.GetHashCode();
        }
        /// <summary>
        /// Converts RangeDouble to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", min, max);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type RangeDouble.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is RangeDouble) ? (this == (RangeDouble)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two RangeDouble objects are equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(RangeDouble a, RangeDouble b)
        {
            return (a.Max == b.Max && a.Min == b.Min);
        }
        /// <summary>
        /// Checks if two RangeDouble objects are not equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(RangeDouble a, RangeDouble b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of RangeDouble.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        object ICloneable.Clone()
        {
            return new RangeDouble(min, max);
        }
        /// <summary>
        /// Creates a copy of RangeDouble.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        public RangeDouble Clone()
        {
            return new RangeDouble(min, max);
        }
        #endregion
    }
}
