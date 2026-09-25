using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of float numbers representing an ordered pair of width and height.
    /// </summary>
    [Serializable]
    public struct SizeFloat : IEquatable<SizeFloat>, ICloneable
    {
        #region Private data
        private float width;
        private float height;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of float numbers representing an ordered pair of width and height.
        /// </summary>
        /// <param name="width">Width.</param>
        /// <param name="height">Height.</param>
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

        #region Size operations
        /// <summary>
        /// Represents a size whose width and height are zero.
        /// </summary>
        public static readonly SizeFloat Empty;

        /// <summary>
        /// Tests whether both dimensions are zero.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public bool IsEmpty => width == 0 && height == 0;

        /// <summary>
        /// Initializes a size from the coordinates of a point.
        /// </summary>
        /// <param name="pt">Point.</param>
        public SizeFloat(PointFloat pt) : this(pt.X, pt.Y)
        {
        }

        /// <summary>
        /// Initializes a size from another size.
        /// </summary>
        /// <param name="size">Size to copy.</param>
        public SizeFloat(SizeFloat size) : this(size.width, size.height)
        {
        }

        /// <summary>
        /// Initializes a size from vector components.
        /// </summary>
        /// <param name="vector">Width and height.</param>
        public SizeFloat(System.Numerics.Vector2 vector) : this(vector.X, vector.Y)
        {
        }

        /// <summary>
        /// Returns the width and height as a vector.
        /// </summary>
        /// <returns>Two-component vector.</returns>
        public System.Numerics.Vector2 ToVector2() => new System.Numerics.Vector2(width, height);

        /// <summary>
        /// Converts a size to a vector.
        /// </summary>
        /// <param name="size">Size.</param>
        public static explicit operator System.Numerics.Vector2(SizeFloat size) => size.ToVector2();

        /// <summary>
        /// Converts a vector to a size.
        /// </summary>
        /// <param name="vector">Width and height.</param>
        public static explicit operator SizeFloat(System.Numerics.Vector2 vector) => new SizeFloat(vector);

        /// <summary>
        /// Returns a point whose coordinates are the width and height.
        /// </summary>
        /// <returns>Point.</returns>
        public PointFloat ToPointF() => (PointFloat)this;

        /// <summary>
        /// Truncates both dimensions toward zero.
        /// </summary>
        /// <returns>Integer size.</returns>
        public SizeInt ToSize() => SizeInt.Truncate(this);

        /// <summary>
        /// Converts the width and height to point coordinates.
        /// </summary>
        /// <param name="size">Size.</param>
        public static explicit operator PointFloat(SizeFloat size) => new PointFloat(size.Width, size.Height);

        /// <summary>
        /// Adds the corresponding dimensions of two sizes.
        /// </summary>
        /// <param name="sz1">First size.</param>
        /// <param name="sz2">Second size.</param>
        /// <returns>Resulting size.</returns>
        public static SizeFloat Add(SizeFloat sz1, SizeFloat sz2) =>
            new SizeFloat(sz1.Width + sz2.Width, sz1.Height + sz2.Height);

        /// <summary>
        /// Adds the corresponding dimensions of two sizes.
        /// </summary>
        /// <param name="sz1">First size.</param>
        /// <param name="sz2">Second size.</param>
        /// <returns>Resulting size.</returns>
        public static SizeFloat operator +(SizeFloat sz1, SizeFloat sz2) => Add(sz1, sz2);

        /// <summary>
        /// Subtracts the corresponding dimensions of two sizes.
        /// </summary>
        /// <param name="sz1">First size.</param>
        /// <param name="sz2">Second size.</param>
        /// <returns>Resulting size.</returns>
        public static SizeFloat Subtract(SizeFloat sz1, SizeFloat sz2) =>
            new SizeFloat(sz1.Width - sz2.Width, sz1.Height - sz2.Height);

        /// <summary>
        /// Subtracts the corresponding dimensions of two sizes.
        /// </summary>
        /// <param name="sz1">First size.</param>
        /// <param name="sz2">Second size.</param>
        /// <returns>Resulting size.</returns>
        public static SizeFloat operator -(SizeFloat sz1, SizeFloat sz2) => Subtract(sz1, sz2);

        /// <summary>
        /// Multiplies both dimensions by a scalar.
        /// </summary>
        /// <param name="left">Multiplier.</param>
        /// <param name="right">Size.</param>
        /// <returns>Scaled size.</returns>
        public static SizeFloat operator *(float left, SizeFloat right) => Multiply(right, left);

        /// <summary>
        /// Multiplies both dimensions by a scalar.
        /// </summary>
        /// <param name="left">Size.</param>
        /// <param name="right">Multiplier.</param>
        /// <returns>Scaled size.</returns>
        public static SizeFloat operator *(SizeFloat left, float right) => Multiply(left, right);

        /// <summary>
        /// Divides both dimensions by a scalar.
        /// </summary>
        /// <param name="left">Size.</param>
        /// <param name="right">Divisor.</param>
        /// <returns>Scaled size.</returns>
        public static SizeFloat operator /(SizeFloat left, float right) =>
            new SizeFloat(left.width / right, left.height / right);

        /// <summary>
        /// Tests whether another size has the same dimensions.
        /// </summary>
        /// <param name="other">Size.</param>
        /// <returns>Boolean.</returns>
        public bool Equals(SizeFloat other) => this == other;

        #endregion

        #region Overrides
        /// <summary>
        /// Returns the hash code for this object.
        /// </summary>
        /// <returns>Integer number.</returns>
        public override int GetHashCode()
        {
            return new System.Drawing.SizeF(width, height).GetHashCode();
        }
        /// <summary>
        /// Converts a SizeFloat to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters.</returns>
        public override string ToString()
        {
            return new System.Drawing.SizeF(width, height).ToString();
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type SizeFloat.
        /// </summary>
        /// <param name="obj">Object.</param>
        /// <returns>Boolean.</returns>
        public override bool Equals(object obj)
        {
            return obj is SizeFloat other && Equals(other);
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two SizeFloat objects are equal.
        /// </summary>
        /// <param name="sz1">Pair of numbers.</param>
        /// <param name="sz2">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator ==(SizeFloat sz1, SizeFloat sz2)
        {
            return (sz1.Width == sz2.Width && sz1.Height == sz2.Height);
        }
        /// <summary>
        /// Checks if two SizeFloat objects are not equal.
        /// </summary>
        /// <param name="sz1">Pair of numbers.</param>
        /// <param name="sz2">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator !=(SizeFloat sz1, SizeFloat sz2)
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
        private static SizeFloat Multiply(SizeFloat size, float multiplier) =>
            new SizeFloat(size.width * multiplier, size.height * multiplier);

        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of SizeFloat.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        object ICloneable.Clone()
        {
            return new SizeFloat(width, height);
        }
        /// <summary>
        /// Creates a copy of SizeFloat.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        public SizeFloat Clone()
        {
            return new SizeFloat(width, height);
        }
        #endregion
    }
}
