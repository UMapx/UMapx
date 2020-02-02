// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;

namespace UMapx.Core
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                                 UMAPX.CORE
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Structures: Range, Point, Size
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
    /// <summary>
    /// Defines a pair of integer numbers representing a line segment.
    /// </summary>
    [Serializable]
    public struct RangeInt : ICloneable
    {
        #region Private data
        private int max;
        private int min;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of integer numbers representing a line segment.
        /// </summary>
        /// <param name="min">Lower bound of the segment</param>
        /// <param name="max">Upper bound of the segment</param>
        public RangeInt(int min, int max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Gets or sets the lower bound of the line segment.
        /// </summary>
        public int Min
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
        public int Max
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
        public bool IsOnRange(int x)
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
        /// Converts RangeInt to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", min, max);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type RangeInt.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is RangeInt) ? (this == (RangeInt)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two RangeDouble objects are equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(RangeInt a, RangeInt b)
        {
            return (a.Max == b.Max && a.Min == b.Min);
        }
        /// <summary>
        /// Checks if two RangeDouble objects are not equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(RangeInt a, RangeInt b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of RangeInt.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        object ICloneable.Clone()
        {
            return new RangeInt(min, max);
        }
        /// <summary>
        /// Creates a copy of RangeInt.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        public RangeInt Clone()
        {
            return new RangeInt(min, max);
        }
        #endregion
    }
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
    /// <summary>
    /// Defines a pair of integer numbers representing an ordered pair of X and Y coordinates.
    /// </summary>
    [Serializable]
    public struct PointInt : ICloneable
    {
        #region Private data
        private int y;
        private int x;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of integer numbers representing an ordered pair of X and Y coordinates.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        public PointInt(int x, int y)
        {
            this.x = x;
            this.y = y;
        }
        /// <summary>
        /// Gets or sets the coordinate X.
        /// </summary>
        public int X
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
        public int Y
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
        /// Converts a PointInt to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", x, y);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type PointInt.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is PointInt) ? (this == (PointInt)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two PointDouble objects are equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(PointInt a, PointInt b)
        {
            return (a.X == b.X && a.Y == b.Y);
        }
        /// <summary>
        /// Checks if two PointDouble objects are not equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(PointInt a, PointInt b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of PointInt.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        object ICloneable.Clone()
        {
            return new PointInt(x, y);
        }
        /// <summary>
        /// Creates a copy of PointInt.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        public PointInt Clone()
        {
            return new PointInt(x, y);
        }
        #endregion
    }
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
    /// <summary>
    /// Defines a pair of integer numbers representing an ordered pair of width and height.
    /// </summary>
    [Serializable]
    public struct SizeInt : ICloneable
    {
        #region Private data
        private int width;
        private int height;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of integer numbers representing an ordered pair of width and height.
        /// </summary>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
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
        /// Converts a SizeInt to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", width, height);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type SizeInt.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is SizeInt) ? (this == (SizeInt)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two SizeInt objects are equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(SizeInt a, SizeInt b)
        {
            return (a.Width == b.Width && a.Height == b.Height);
        }
        /// <summary>
        /// Checks if two SizeInt objects are not equal.
        /// </summary>
        /// <param name="a">Pair of numbers</param>
        /// <param name="b">Pair of numbers</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(SizeInt a, SizeInt b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of SizeInt.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        object ICloneable.Clone()
        {
            return new SizeInt(width, height);
        }
        /// <summary>
        /// Creates a copy of SizeInt.
        /// </summary>
        /// <returns>Pair of numbers</returns>
        public SizeInt Clone()
        {
            return new SizeInt(width, height);
        }
        #endregion
    }
    #endregion
}
