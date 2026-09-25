using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of integer numbers representing an ordered pair of width and height.
    /// </summary>
    [Serializable]
    public struct SizeInt : IEquatable<SizeInt>, ICloneable
    {
        #region Private data
        private int width;
        private int height;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of integer numbers representing an ordered pair of width and height.
        /// </summary>
        /// <param name="width">Width.</param>
        /// <param name="height">Height.</param>
        public SizeInt(int width, int height)
        {
            this.height = height;
            this.width = width;
        }
        /// <summary>
        /// Gets or sets the height.
        /// </summary>
        public int Height
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
        public int Width
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

        #region Size operations
        /// <summary>
        /// Represents a size whose width and height are zero.
        /// </summary>
        public static readonly SizeInt Empty;

        /// <summary>
        /// Tests whether both dimensions are zero.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public bool IsEmpty => width == 0 && height == 0;

        /// <summary>
        /// Initializes a size from the coordinates of a point.
        /// </summary>
        /// <param name="pt">Point.</param>
        public SizeInt(PointInt pt) : this(pt.X, pt.Y)
        {
        }

        /// <summary>
        /// Converts integer dimensions to single precision.
        /// </summary>
        /// <param name="p">Integer size.</param>
        public static implicit operator SizeFloat(SizeInt p) => new SizeFloat(p.Width, p.Height);

        /// <summary>
        /// Rounds both dimensions toward positive infinity.
        /// </summary>
        /// <param name="value">Floating-point size.</param>
        /// <returns>Integer size.</returns>
        public static SizeInt Ceiling(SizeFloat value) =>
            new SizeInt(unchecked((int)Math.Ceiling(value.Width)), unchecked((int)Math.Ceiling(value.Height)));

        /// <summary>
        /// Rounds both dimensions to the nearest integer, with midpoint ties to even.
        /// </summary>
        /// <param name="value">Floating-point size.</param>
        /// <returns>Integer size.</returns>
        public static SizeInt Round(SizeFloat value) =>
            new SizeInt(unchecked((int)Math.Round(value.Width)), unchecked((int)Math.Round(value.Height)));

        /// <summary>
        /// Truncates both dimensions toward zero.
        /// </summary>
        /// <param name="value">Floating-point size.</param>
        /// <returns>Integer size.</returns>
        public static SizeInt Truncate(SizeFloat value) =>
            new SizeInt(unchecked((int)value.Width), unchecked((int)value.Height));

        /// <summary>
        /// Converts the width and height to point coordinates.
        /// </summary>
        /// <param name="size">Size.</param>
        public static explicit operator PointInt(SizeInt size) => new PointInt(size.Width, size.Height);

        /// <summary>
        /// Adds the corresponding dimensions of two sizes.
        /// </summary>
        /// <param name="sz1">First size.</param>
        /// <param name="sz2">Second size.</param>
        /// <returns>Resulting size.</returns>
        public static SizeInt Add(SizeInt sz1, SizeInt sz2) =>
            new SizeInt(unchecked(sz1.Width + sz2.Width), unchecked(sz1.Height + sz2.Height));

        /// <summary>
        /// Adds the corresponding dimensions of two sizes.
        /// </summary>
        /// <param name="sz1">First size.</param>
        /// <param name="sz2">Second size.</param>
        /// <returns>Resulting size.</returns>
        public static SizeInt operator +(SizeInt sz1, SizeInt sz2) => Add(sz1, sz2);

        /// <summary>
        /// Subtracts the corresponding dimensions of two sizes.
        /// </summary>
        /// <param name="sz1">First size.</param>
        /// <param name="sz2">Second size.</param>
        /// <returns>Resulting size.</returns>
        public static SizeInt Subtract(SizeInt sz1, SizeInt sz2) =>
            new SizeInt(unchecked(sz1.Width - sz2.Width), unchecked(sz1.Height - sz2.Height));

        /// <summary>
        /// Subtracts the corresponding dimensions of two sizes.
        /// </summary>
        /// <param name="sz1">First size.</param>
        /// <param name="sz2">Second size.</param>
        /// <returns>Resulting size.</returns>
        public static SizeInt operator -(SizeInt sz1, SizeInt sz2) => Subtract(sz1, sz2);

        /// <summary>
        /// Multiplies both dimensions by a scalar.
        /// </summary>
        /// <param name="left">Multiplier.</param>
        /// <param name="right">Size.</param>
        /// <returns>Scaled size.</returns>
        public static SizeInt operator *(int left, SizeInt right) => Multiply(right, left);

        /// <summary>
        /// Multiplies both dimensions by a scalar.
        /// </summary>
        /// <param name="left">Size.</param>
        /// <param name="right">Multiplier.</param>
        /// <returns>Scaled size.</returns>
        public static SizeInt operator *(SizeInt left, int right) => Multiply(left, right);

        /// <summary>
        /// Divides both dimensions by a scalar.
        /// </summary>
        /// <param name="left">Size.</param>
        /// <param name="right">Divisor.</param>
        /// <returns>Scaled size.</returns>
        public static SizeInt operator /(SizeInt left, int right) =>
            new SizeInt(unchecked(left.width / right), unchecked(left.height / right));

        /// <summary>
        /// Multiplies both dimensions by a scalar.
        /// </summary>
        /// <param name="left">Multiplier.</param>
        /// <param name="right">Size.</param>
        /// <returns>Scaled size.</returns>
        public static SizeFloat operator *(float left, SizeInt right) => Multiply(right, left);

        /// <summary>
        /// Multiplies both dimensions by a scalar.
        /// </summary>
        /// <param name="left">Size.</param>
        /// <param name="right">Multiplier.</param>
        /// <returns>Scaled size.</returns>
        public static SizeFloat operator *(SizeInt left, float right) => Multiply(left, right);

        /// <summary>
        /// Divides both dimensions by a scalar.
        /// </summary>
        /// <param name="left">Size.</param>
        /// <param name="right">Divisor.</param>
        /// <returns>Scaled size.</returns>
        public static SizeFloat operator /(SizeInt left, float right) =>
            new SizeFloat(left.width / right, left.height / right);

        /// <summary>
        /// Tests whether another size has the same dimensions.
        /// </summary>
        /// <param name="other">Size.</param>
        /// <returns>Boolean.</returns>
        public bool Equals(SizeInt other) => this == other;

        #endregion

        #region Overrides
        /// <summary>
        /// Returns the hash code for this object.
        /// </summary>
        /// <returns>Integer number.</returns>
        public override int GetHashCode()
        {
            return unchecked((width * 397) ^ height);
        }
        /// <summary>
        /// Converts a SizeInt to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters.</returns>
        public override string ToString()
        {
            return string.Format("{{Width={0}, Height={1}}}", width, height);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type SizeInt.
        /// </summary>
        /// <param name="obj">Object.</param>
        /// <returns>Boolean.</returns>
        public override bool Equals(object obj)
        {
            return obj is SizeInt other && Equals(other);
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two SizeInt objects are equal.
        /// </summary>
        /// <param name="sz1">Pair of numbers.</param>
        /// <param name="sz2">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator ==(SizeInt sz1, SizeInt sz2)
        {
            return (sz1.Width == sz2.Width && sz1.Height == sz2.Height);
        }
        /// <summary>
        /// Checks if two SizeInt objects are not equal.
        /// </summary>
        /// <param name="sz1">Pair of numbers.</param>
        /// <param name="sz2">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator !=(SizeInt sz1, SizeInt sz2)
        {
            return !(sz1 == sz2);
        }
        #endregion

        #region Private methods
        /// <summary>
        /// Multiplies both dimensions by a scalar.
        /// </summary>
        /// <param name="size">Size.</param>
        /// <param name="multiplier">Multiplier.</param>
        /// <returns>Scaled size.</returns>
        private static SizeInt Multiply(SizeInt size, int multiplier) =>
            new SizeInt(unchecked(size.width * multiplier), unchecked(size.height * multiplier));

        /// <summary>
        /// Multiplies both dimensions by a scalar.
        /// </summary>
        /// <param name="size">Size.</param>
        /// <param name="multiplier">Multiplier.</param>
        /// <returns>Scaled size.</returns>
        private static SizeFloat Multiply(SizeInt size, float multiplier) =>
            new SizeFloat(size.width * multiplier, size.height * multiplier);

        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of SizeInt.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        object ICloneable.Clone()
        {
            return new SizeInt(width, height);
        }
        /// <summary>
        /// Creates a copy of SizeInt.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        public SizeInt Clone()
        {
            return new SizeInt(width, height);
        }
        #endregion
    }
}
