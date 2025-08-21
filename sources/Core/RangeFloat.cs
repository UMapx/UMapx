using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of float numbers representing a line segment.
    /// </summary>
    [Serializable]
    public struct RangeFloat : ICloneable
    {
        #region Private data
        private float max;
        private float min;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of float numbers representing a line segment.
        /// </summary>
        /// <param name="min">Lower bound of the segment</param>
        /// <param name="max">Upper bound of the segment</param>
        public RangeFloat(float min, float max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Gets or sets the lower bound of the line segment.
        /// </summary>
        public float Min
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
        public float Max
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
        /// <param name="x">Value</param>
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
        /// Converts RangeFloat to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", min, max);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type RangeFloat.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is RangeFloat) ? (this == (RangeFloat)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two RangeFloat objects are equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(RangeFloat a, RangeFloat b)
        {
            return (a.Max == b.Max && a.Min == b.Min);
        }
        /// <summary>
        /// Checks if two RangeFloat objects are not equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(RangeFloat a, RangeFloat b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of RangeFloat.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        object ICloneable.Clone()
        {
            return new RangeFloat(min, max);
        }
        /// <summary>
        /// Creates a copy of RangeFloat.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        public RangeFloat Clone()
        {
            return new RangeFloat(min, max);
        }
        #endregion
    }
}
