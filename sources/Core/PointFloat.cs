using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a pair of float numbers representing an ordered pair of X and Y coordinates.
    /// </summary>
    [Serializable]
    public struct PointFloat : IEquatable<PointFloat>, ICloneable
    {
        #region Private data
        private float y;
        private float x;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes a pair of float numbers representing an ordered pair of X and Y coordinates.
        /// </summary>
        /// <param name="x">Coordinate X.</param>
        /// <param name="y">Coordinate Y.</param>
        public PointFloat(float x, float y)
        {
            this.x = x;
            this.y = y;
        }
        /// <summary>
        /// Gets or sets the coordinate X.
        /// </summary>
        public float X
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
        public float Y
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
        public static readonly PointFloat Empty;

        /// <summary>
        /// Tests whether both coordinates are zero.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public bool IsEmpty => x == 0 && y == 0;

        /// <summary>
        /// Initializes a point from a two-component vector.
        /// </summary>
        /// <param name="vector">Coordinates.</param>
        public PointFloat(System.Numerics.Vector2 vector) : this(vector.X, vector.Y)
        {
        }

        /// <summary>
        /// Returns a vector containing the X and Y coordinates.
        /// </summary>
        /// <returns>Coordinate vector.</returns>
        public System.Numerics.Vector2 ToVector2() => new System.Numerics.Vector2(x, y);

        /// <summary>
        /// Converts a point to a coordinate vector.
        /// </summary>
        /// <param name="point">Point.</param>
        public static explicit operator System.Numerics.Vector2(PointFloat point) => point.ToVector2();

        /// <summary>
        /// Converts a coordinate vector to a point.
        /// </summary>
        /// <param name="vector">Coordinates.</param>
        public static explicit operator PointFloat(System.Numerics.Vector2 vector) => new PointFloat(vector);

        /// <summary>
        /// Returns a point translated by the specified dimensions.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Horizontal and vertical displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointFloat Add(PointFloat pt, SizeInt sz)
        {
            return new PointFloat(pt.x + sz.Width, pt.y + sz.Height);
        }

        /// <summary>
        /// Returns a point translated by the negative of the specified dimensions.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Horizontal and vertical displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointFloat Subtract(PointFloat pt, SizeInt sz)
        {
            return new PointFloat(pt.x - sz.Width, pt.y - sz.Height);
        }

        /// <summary>
        /// Adds a displacement to a point.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointFloat operator +(PointFloat pt, SizeInt sz) => Add(pt, sz);

        /// <summary>
        /// Subtracts a displacement from a point.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointFloat operator -(PointFloat pt, SizeInt sz) => Subtract(pt, sz);

        /// <summary>
        /// Returns a point translated by the specified dimensions.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Horizontal and vertical displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointFloat Add(PointFloat pt, SizeFloat sz)
        {
            return new PointFloat(pt.x + sz.Width, pt.y + sz.Height);
        }

        /// <summary>
        /// Returns a point translated by the negative of the specified dimensions.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Horizontal and vertical displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointFloat Subtract(PointFloat pt, SizeFloat sz)
        {
            return new PointFloat(pt.x - sz.Width, pt.y - sz.Height);
        }

        /// <summary>
        /// Adds a displacement to a point.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointFloat operator +(PointFloat pt, SizeFloat sz) => Add(pt, sz);

        /// <summary>
        /// Subtracts a displacement from a point.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <param name="sz">Displacement.</param>
        /// <returns>Translated point.</returns>
        public static PointFloat operator -(PointFloat pt, SizeFloat sz) => Subtract(pt, sz);

        /// <summary>
        /// Tests whether another point has the same coordinates.
        /// </summary>
        /// <param name="other">Point.</param>
        /// <returns>Boolean.</returns>
        public bool Equals(PointFloat other) => this == other;
        #endregion

        #region Point arithmetic
        /// <summary>
        /// Translates each point by adding the specified coordinate offsets.
        /// </summary>
        /// <param name="points">Points to translate.</param>
        /// <param name="point">Horizontal and vertical offsets to add.</param>
        /// <returns>New array containing the translated points in the original order.</returns>
        public static PointFloat[] Add(PointFloat[] points, PointFloat point)
        {
            var count = points.Length;
            var output = new PointFloat[count];

            for (int i = 0; i < count; i++)
            {
                output[i] = new PointFloat
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
        public static PointFloat[] Sub(PointFloat[] points, PointFloat point)
        {
            var count = points.Length;
            var output = new PointFloat[count];

            for (int i = 0; i < count; i++)
            {
                output[i] = new PointFloat
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
        public static PointFloat[] Rotate(PointFloat[] points, PointFloat centerPoint, float angle)
        {
            int length = points.Length;
            var output = new PointFloat[length];

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
        /// <returns>Rotated point with single-precision floating-point coordinates.</returns>
        public PointFloat Rotate(PointFloat centerPoint, float angleInDegrees)
        {
            var pointToRotate = this;
            double angleInRadians = angleInDegrees * (Math.PI / 180);
            double cosTheta = Math.Cos(angleInRadians);
            double sinTheta = Math.Sin(angleInRadians);

            return new PointFloat
            {
                X =
                    (float)
                    (cosTheta * (pointToRotate.X - centerPoint.X) -
                    sinTheta * (pointToRotate.Y - centerPoint.Y) + centerPoint.X),
                Y =
                    (float)
                    (sinTheta * (pointToRotate.X - centerPoint.X) +
                    cosTheta * (pointToRotate.Y - centerPoint.Y) + centerPoint.Y)
            };
        }

        /// <summary>
        /// Calculates the axis-aligned coordinate bounds of a point array.
        /// </summary>
        /// <param name="points">Nonempty array of points with finite coordinates to bound.</param>
        /// <returns>Rectangle whose left and top edges are the minimum coordinates and whose width and height are the coordinate spans.</returns>
        public static RectangleFloat GetRectangle(PointFloat[] points)
        {
            int length = points.Length;
            float xmin = float.MaxValue;
            float ymin = float.MaxValue;
            float xmax = float.MinValue;
            float ymax = float.MinValue;

            for (int i = 0; i < length; i++)
            {
                float x = points[i].X;
                float y = points[i].Y;

                if (x < xmin)
                    xmin = x;
                if (y < ymin)
                    ymin = y;
                if (x > xmax)
                    xmax = x;
                if (y > ymax)
                    ymax = y;
            }

            return new RectangleFloat(xmin, ymin, xmax - xmin, ymax - ymin);
        }

        /// <summary>
        /// Calculates a signed angle at the left point using the support and right points.
        /// </summary>
        /// <param name="right">Point defining one ray from the vertex.</param>
        /// <param name="support">Point defining the other ray from the vertex.</param>
        /// <returns>Signed angle in approximate degrees, or NaN if the argument passed to the inverse cosine is outside its valid range.</returns>
        public float GetAngle(PointFloat right, PointFloat support)
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
        /// <param name="right">Point supplying the X coordinate.</param>
        /// <returns>Point with coordinates <c>(right.X, left.Y)</c>.</returns>
        public PointFloat GetSupportedPoint(PointFloat right)
        {
            var left = this;
            return new PointFloat(right.X, left.Y);
        }

        /// <summary>
        /// Calculates the arithmetic mean of the point coordinates.
        /// </summary>
        /// <param name="points">Points to average.</param>
        /// <returns>Point containing the mean X and Y coordinates, or NaN in both coordinates for an empty array.</returns>
        public static PointFloat GetMeanPoint(params PointFloat[] points)
        {
            var point = new PointFloat(0, 0);
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
            return new System.Drawing.PointF(x, y).GetHashCode();
        }
        /// <summary>
        /// Converts a PointFloat to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters.</returns>
        public override string ToString()
        {
            return new System.Drawing.PointF(x, y).ToString();
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type PointFloat.
        /// </summary>
        /// <param name="obj">Object.</param>
        /// <returns>Boolean.</returns>
        public override bool Equals(object obj)
        {
            return obj is PointFloat other && Equals(other);
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two PointFloat objects are equal.
        /// </summary>
        /// <param name="left">Pair of numbers.</param>
        /// <param name="right">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator ==(PointFloat left, PointFloat right)
        {
            return (left.X == right.X && left.Y == right.Y);
        }
        /// <summary>
        /// Checks if two PointFloat objects are not equal.
        /// </summary>
        /// <param name="left">Pair of numbers.</param>
        /// <param name="right">Pair of numbers.</param>
        /// <returns>Boolean.</returns>
        public static bool operator !=(PointFloat left, PointFloat right)
        {
            return !(left == right);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of PointFloat.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        object ICloneable.Clone()
        {
            return new PointFloat(x, y);
        }
        /// <summary>
        /// Creates a copy of PointFloat.
        /// </summary>
        /// <returns>Pair of numbers.</returns>
        public PointFloat Clone()
        {
            return new PointFloat(x, y);
        }
        #endregion
    }
}
