using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of integer numbers representing an ordered pair of X and Y coordinates.
    /// </summary>
    [Serializable]
    public struct PointInt : IEquatable<PointInt>, ICloneable
    {
        #region Private data
        private int y;
        private int x;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of integer numbers representing an ordered pair of X and Y coordinates.
        /// </summary>
        /// <param name="x">Coordinate X.</param>
        /// <param name="y">Coordinate Y.</param>
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

        #region Standard point operations
        /// <summary>
        /// Represents the origin, with both coordinates equal to zero.
        /// </summary>
        public static readonly PointInt Empty;

        /// <summary>
        /// Tests whether both coordinates are zero.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public bool IsEmpty => x == 0 && y == 0;

        /// <summary>
        /// Initializes a point whose coordinates are the specified width and height.
        /// </summary>
        /// <param name="sz">Coordinates represented as dimensions.</param>
        public PointInt(SizeInt sz) : this(sz.Width, sz.Height)
        {
        }

        /// <summary>
        /// Initializes a point from two packed signed 16-bit coordinates.
        /// </summary>
        /// <param name="dw">Low 16 bits contain X; high 16 bits contain Y.</param>
        public PointInt(int dw) : this(unchecked((short)dw), unchecked((short)(dw >> 16)))
        {
        }

        /// <summary>
        /// Converts integer coordinates to single precision.
        /// </summary>
        /// <param name="p">Integer point.</param>
        public static implicit operator PointFloat(PointInt p) => new PointFloat(p.x, p.y);

        /// <summary>
        /// Converts coordinates to width and height.
        /// </summary>
        /// <param name="p">Point.</param>
        public static explicit operator SizeInt(PointInt p) => new SizeInt(p.x, p.y);

        /// <summary>
        /// Rounds each coordinate toward positive infinity.
        /// </summary>
        /// <param name="value">Floating-point point.</param>
        /// <returns>Integer point.</returns>
        public static PointInt Ceiling(PointFloat value)
        {
            return new PointInt(unchecked((int)Math.Ceiling(value.X)), unchecked((int)Math.Ceiling(value.Y)));
        }

        /// <summary>
        /// Rounds each coordinate to the nearest integer, with midpoint ties to even.
        /// </summary>
        /// <param name="value">Floating-point point.</param>
        /// <returns>Integer point.</returns>
        public static PointInt Round(PointFloat value)
        {
            return new PointInt(unchecked((int)Math.Round(value.X)), unchecked((int)Math.Round(value.Y)));
        }

        /// <summary>
        /// Truncates each coordinate toward zero.
        /// </summary>
        /// <param name="value">Floating-point point.</param>
        /// <returns>Integer point.</returns>
        public static PointInt Truncate(PointFloat value)
        {
            return new PointInt(unchecked((int)value.X), unchecked((int)value.Y));
        }

        /// <summary>
        /// Returns a point translated by the specified dimensions.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Horizontal and vertical displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointInt Add(PointInt pt, SizeInt sz)
        {
            return new PointInt(unchecked(pt.x + sz.Width), unchecked(pt.y + sz.Height));
        }

        /// <summary>
        /// Returns a point translated by the negative of the specified dimensions.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Horizontal and vertical displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointInt Subtract(PointInt pt, SizeInt sz)
        {
            return new PointInt(unchecked(pt.x - sz.Width), unchecked(pt.y - sz.Height));
        }

        /// <summary>
        /// Adds a displacement to a point.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointInt operator +(PointInt pt, SizeInt sz) => Add(pt, sz);

        /// <summary>
        /// Subtracts a displacement from a point.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointInt operator -(PointInt pt, SizeInt sz) => Subtract(pt, sz);

        /// <summary>
        /// Moves this point by the specified displacement.
        /// </summary>
        /// <param name="dx">Horizontal displacement.</param>
        /// <param name="dy">Vertical displacement.</param>
        public void Offset(int dx, int dy)
        {
            x = unchecked(x + dx);
            y = unchecked(y + dy);
        }

        /// <summary>
        /// Moves this point by the specified displacement.
        /// </summary>
        /// <param name="p">Displacement.</param>
        public void Offset(PointInt p) => Offset(p.x, p.y);

        /// <summary>
        /// Tests whether another point has the same coordinates.
        /// </summary>
        /// <param name="other">Point.</param>
        /// <returns>Boolean.</returns>
        public bool Equals(PointInt other) => this == other;
        #endregion

        #region Point arithmetic
        /// <summary>
        /// Translates each point by adding the specified coordinate offsets.
        /// </summary>
        /// <param name="points">Points to translate.</param>
        /// <param name="point">Horizontal and vertical offsets to add.</param>
        /// <returns>New array containing the translated points in the original order.</returns>
        public static PointInt[] Add(PointInt[] points, PointInt point)
        {
            var count = points.Length;
            var output = new PointInt[count];

            for (int i = 0; i < count; i++)
            {
                output[i] = new PointInt
                {
                    X = points[i].X + point.X,
                    Y = points[i].Y + point.Y
                };
            }

            return output;
        }

        /// <summary>
        /// Translates each point by subtracting the specified coordinate offsets.
        /// </summary>
        /// <param name="points">Points to translate.</param>
        /// <param name="point">Horizontal and vertical offsets to subtract.</param>
        /// <returns>New array containing the translated points in the original order.</returns>
        public static PointInt[] Sub(PointInt[] points, PointInt point)
        {
            var count = points.Length;
            var output = new PointInt[count];

            for (int i = 0; i < count; i++)
            {
                output[i] = new PointInt
                {
                    X = points[i].X - point.X,
                    Y = points[i].Y - point.Y
                };
            }

            return output;
        }

        /// <summary>
        /// Rotates each point around the specified center.
        /// </summary>
        /// <param name="points">Points to rotate.</param>
        /// <param name="centerPoint">Center of rotation.</param>
        /// <param name="angle">Rotation angle in degrees.</param>
        /// <returns>New array containing the rotated points in the original order.</returns>
        public static PointInt[] Rotate(PointInt[] points, PointInt centerPoint, float angle)
        {
            int length = points.Length;
            var output = new PointInt[length];

            for (int i = 0; i < length; i++)
            {
                output[i] = points[i].Rotate(centerPoint, angle);
            }

            return output;
        }

        /// <summary>
        /// Rotates a point around the specified center.
        /// </summary>
        /// <param name="centerPoint">Center of rotation.</param>
        /// <param name="angleInDegrees">Rotation angle in degrees.</param>
        /// <returns>Rotated point with coordinates truncated toward zero to integers.</returns>
        public PointInt Rotate(PointInt centerPoint, float angleInDegrees)
        {
            var pointToRotate = this;
            double angleInRadians = angleInDegrees * (Math.PI / 180);
            double cosTheta = Math.Cos(angleInRadians);
            double sinTheta = Math.Sin(angleInRadians);

            return new PointInt
            {
                X =
                    (int)
                    (cosTheta * (pointToRotate.X - centerPoint.X) -
                    sinTheta * (pointToRotate.Y - centerPoint.Y) + centerPoint.X),
                Y =
                    (int)
                    (sinTheta * (pointToRotate.X - centerPoint.X) +
                    cosTheta * (pointToRotate.Y - centerPoint.Y) + centerPoint.Y)
            };
        }

        /// <summary>
        /// Calculates the axis-aligned coordinate bounds of a point array.
        /// </summary>
        /// <param name="points">Nonempty array of points to bound.</param>
        /// <returns>RectangleInt whose left and top edges are the minimum coordinates and whose width and height are the coordinate spans.</returns>
        public static RectangleInt GetRectangle(PointInt[] points)
        {
            int length = points.Length;
            int xmin = int.MaxValue;
            int ymin = int.MaxValue;
            int xmax = int.MinValue;
            int ymax = int.MinValue;

            for (int i = 0; i < length; i++)
            {
                int x = points[i].X;
                int y = points[i].Y;

                if (x < xmin)
                    xmin = x;
                if (y < ymin)
                    ymin = y;
                if (x > xmax)
                    xmax = x;
                if (y > ymax)
                    ymax = y;
            }

            return new RectangleInt(xmin, ymin, xmax - xmin, ymax - ymin);
        }

        /// <summary>
        /// Calculates a signed angle at the left point using the support and right points.
        /// </summary>
        /// <param name="right">PointInt defining one ray from the vertex.</param>
        /// <param name="support">PointInt defining the other ray from the vertex.</param>
        /// <returns>Signed angle in approximate degrees, or NaN if the argument passed to the inverse cosine is outside its valid range.</returns>
        public float GetAngle(PointInt right, PointInt support)
        {
            var left = this;
            double kk = left.Y > right.Y ? 1 : -1;

            double x1 = left.X - support.X;
            double y1 = left.Y - support.Y;

            double x2 = right.X - left.X;
            double y2 = right.Y - left.Y;

            double a = Math.Sqrt(x1 * x1 + y1 * y1);
            double b = Math.Sqrt(x2 * x2 + y2 * y2);
            double c = x1 * x2 + y1 * y2;

            double d = Div(Div(c, a), b);

            return (float)(kk * (180.0 - Math.Acos(d) * 57.3));
        }

        /// <summary>
        /// Creates a point by combining coordinates from two points.
        /// </summary>
        /// <param name="right">PointInt supplying the X coordinate.</param>
        /// <returns>PointInt with coordinates <c>(right.X, left.Y)</c>.</returns>
        public PointInt GetSupportedPoint(PointInt right)
        {
            var left = this;
            return new PointInt(right.X, left.Y);
        }

        /// <summary>
        /// Calculates the arithmetic mean of the point coordinates.
        /// </summary>
        /// <param name="points">Nonempty array of points to average.</param>
        /// <returns>PointInt containing the mean X and Y coordinates, with integer division truncating each result toward zero.</returns>
        /// <exception cref="DivideByZeroException">The array contains no points.</exception>
        public static PointInt GetMeanPoint(params PointInt[] points)
        {
            var point = new PointInt(0, 0);
            var length = points.Length;

            for (int i = 0; i < length; i++)
            {
                point.X += points[i].X;
                point.Y += points[i].Y;
            }

            point.X /= length;
            point.Y /= length;

            return point;
        }
        #endregion

        #region Private methods

        /// <summary>
        /// Divides two values, substituting the smallest positive double value for zero divided by zero.
        /// </summary>
        /// <param name="a">Numerator.</param>
        /// <param name="b">Denominator.</param>
        /// <returns><see cref="double.Epsilon"/> when both operands are zero; otherwise, the result of dividing the numerator by the denominator.</returns>
        private static double Div(double a, double b)
        {
            if (a == 0 && b == 0)
            {
                return double.Epsilon;
            }

            return a / b;
        }

        #endregion

        #region Overrides
        /// <summary>
        /// Returns the hash code for this object.
        /// </summary>
        /// <returns>Integer number.</returns>
        public override int GetHashCode()
        {
            return unchecked((x * 397) ^ y);
        }
        /// <summary>
        /// Converts a PointInt to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters.</returns>
        public override string ToString()
        {
            return string.Format("{{X={0},Y={1}}}", x, y);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type PointInt.
        /// </summary>
        /// <param name="obj">Object.</param>
        /// <returns>Boolean.</returns>
        public override bool Equals(object obj)
        {
            return obj is PointInt other && Equals(other);
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two PointInt objects are equal.
        /// </summary>
        /// <param name="left">Pair of numbers.</param>
        /// <param name="right">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator ==(PointInt left, PointInt right)
        {
            return (left.X == right.X && left.Y == right.Y);
        }
        /// <summary>
        /// Checks if two PointInt objects are not equal.
        /// </summary>
        /// <param name="left">Pair of numbers.</param>
        /// <param name="right">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator !=(PointInt left, PointInt right)
        {
            return !(left == right);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of PointInt.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        object ICloneable.Clone()
        {
            return new PointInt(x, y);
        }
        /// <summary>
        /// Creates a copy of PointInt.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        public PointInt Clone()
        {
            return new PointInt(x, y);
        }
        #endregion
    }
}
