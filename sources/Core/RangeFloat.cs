using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of float numbers representing a line segment.
    /// </summary>
    [Serializable]
    public struct RangeFloat : IEquatable<RangeFloat>, ICloneable
    {
        #region Private data
        private float max;
        private float min;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of float numbers representing a line segment.
        /// </summary>
        /// <param name="min">Lower bound of the segment.</param>
        /// <param name="max">Upper bound of the segment.</param>
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
        /// <param name="x">Value.</param>
        /// <returns>Boolean.</returns>
        public bool IsOnRange(float x)
        {
            if ((x >= this.min) && (x <= this.max))
            {
                return true;
            }
            return false;
        }
        #endregion

        #region Range operations
        /// <summary>
        /// Represents a range with no values.
        /// </summary>
        public static readonly RangeFloat Empty = new RangeFloat(0, -1);

        /// <summary>
        /// Gets the upper bound minus the lower bound.
        /// </summary>
        public float Length => max - min;

        /// <summary>
        /// Tests whether a value lies between the inclusive bounds.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Boolean.</returns>
        public bool Contains(float x) => IsOnRange(x);

        /// <summary>
        /// Tests whether a range lies entirely between the inclusive bounds.
        /// </summary>
        /// <param name="range">Range.</param>
        /// <returns>Boolean.</returns>
        public bool Contains(RangeFloat range)
        {
            return min <= max && range.min <= range.max && min <= range.min && range.max <= max;
        }

        /// <summary>
        /// Tests whether two ranges share at least one value.
        /// </summary>
        /// <param name="range">Range.</param>
        /// <returns>Boolean.</returns>
        public bool IntersectsWith(RangeFloat range)
        {
            return min <= max && range.min <= range.max && min <= range.max && range.min <= max;
        }

        /// <summary>
        /// Replaces this range with its intersection with another range.
        /// </summary>
        /// <param name="range">Range.</param>
        public void Intersect(RangeFloat range)
        {
            this = Intersect(this, range);
        }

        /// <summary>
        /// Returns the intersection of two ranges.
        /// </summary>
        /// <param name="a">First range.</param>
        /// <param name="b">Second range.</param>
        /// <returns>Common range, or Empty if the ranges do not intersect.</returns>
        public static RangeFloat Intersect(RangeFloat a, RangeFloat b)
        {
            if (!a.IntersectsWith(b))
                return Empty;

            return new RangeFloat(Math.Max(a.min, b.min), Math.Min(a.max, b.max));
        }

        /// <summary>
        /// Returns the smallest range containing both ranges.
        /// </summary>
        /// <param name="a">First range.</param>
        /// <param name="b">Second range.</param>
        /// <returns>Bounding range. Ranges with unordered bounds are ignored.</returns>
        public static RangeFloat Union(RangeFloat a, RangeFloat b)
        {
            if (!(a.min <= a.max))
                return b.min <= b.max ? b : Empty;
            if (!(b.min <= b.max))
                return a;

            return new RangeFloat(Math.Min(a.min, b.min), Math.Max(a.max, b.max));
        }

        /// <summary>
        /// Moves both bounds by the specified amount.
        /// </summary>
        /// <param name="offset">Displacement.</param>
        public void Offset(float offset)
        {
            min = min + offset;
            max = max + offset;
        }

        /// <summary>
        /// Expands both ends by the specified amount.
        /// </summary>
        /// <param name="amount">Amount added at each end.</param>
        public void Inflate(float amount)
        {
            min = min - amount;
            max = max + amount;
        }

        /// <summary>
        /// Returns a range expanded at both ends by the specified amount.
        /// </summary>
        /// <param name="range">Range.</param>
        /// <param name="amount">Amount added at each end.</param>
        /// <returns>Expanded range.</returns>
        public static RangeFloat Inflate(RangeFloat range, float amount)
        {
            range.Inflate(amount);
            return range;
        }

        /// <summary>
        /// Tests whether another range has the same bounds.
        /// </summary>
        /// <param name="other">Range.</param>
        /// <returns>Boolean.</returns>
        public bool Equals(RangeFloat other) => this == other;

        #endregion

        #region Overrides
        /// <summary>
        /// Returns the hash code for this object.
        /// </summary>
        /// <returns>Integer number.</returns>
        public override int GetHashCode()
        {
            return min.GetHashCode() ^ max.GetHashCode();
        }
        /// <summary>
        /// Converts RangeFloat to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters.</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", min, max);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type RangeFloat.
        /// </summary>
        /// <param name="obj">Object.</param>
        /// <returns>Boolean.</returns>
        public override bool Equals(object obj)
        {
            return obj is RangeFloat other && Equals(other);
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two RangeFloat objects are equal.
        /// </summary>
        /// <param name="a">Pair of numbers.</param>
        /// <param name="b">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator ==(RangeFloat a, RangeFloat b)
        {
            return (a.Max == b.Max && a.Min == b.Min);
        }
        /// <summary>
        /// Checks if two RangeFloat objects are not equal.
        /// </summary>
        /// <param name="a">Pair of numbers.</param>
        /// <param name="b">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator !=(RangeFloat a, RangeFloat b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of RangeFloat.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        object ICloneable.Clone()
        {
            return new RangeFloat(min, max);
        }
        /// <summary>
        /// Creates a copy of RangeFloat.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        public RangeFloat Clone()
        {
            return new RangeFloat(min, max);
        }
        #endregion
    }
}
