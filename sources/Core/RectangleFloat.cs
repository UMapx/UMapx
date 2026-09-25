using System;
using System.Numerics;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a rectangle with single-precision floating-point coordinates and dimensions.
    /// </summary>
    [Serializable]
    public struct RectangleFloat : IEquatable<RectangleFloat>, ICloneable
    {
        #region Private data
        private float x;
        private float y;
        private float width;
        private float height;
        #endregion

        #region Structure components
        /// <summary>
        /// Represents a rectangle whose coordinates and dimensions are zero.
        /// </summary>
        public static readonly RectangleFloat Empty;

        /// <summary>
        /// Initializes a rectangle with the specified location and size.
        /// </summary>
        /// <param name="x">Left coordinate.</param>
        /// <param name="y">Top coordinate.</param>
        /// <param name="width">Width.</param>
        /// <param name="height">Height.</param>
        public RectangleFloat(float x, float y, float width, float height)
        {
            this.x = x;
            this.y = y;
            this.width = width;
            this.height = height;
        }

        /// <summary>
        /// Initializes a rectangle from a point and a size.
        /// </summary>
        /// <param name="location">Upper-left corner.</param>
        /// <param name="size">Dimensions.</param>
        public RectangleFloat(PointFloat location, SizeFloat size)
            : this(location.X, location.Y, size.Width, size.Height)
        {
        }

        /// <summary>
        /// Gets or sets the left coordinate.
        /// </summary>
        public float X { get => x; set => x = value; }

        /// <summary>
        /// Gets or sets the top coordinate.
        /// </summary>
        public float Y { get => y; set => y = value; }

        /// <summary>
        /// Gets or sets the width.
        /// </summary>
        public float Width { get => width; set => width = value; }

        /// <summary>
        /// Gets or sets the height.
        /// </summary>
        public float Height { get => height; set => height = value; }

        /// <summary>
        /// Gets or sets the upper-left corner.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public PointFloat Location
        {
            get => new PointFloat(x, y);
            set { x = value.X; y = value.Y; }
        }

        /// <summary>
        /// Gets or sets the dimensions.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public SizeFloat Size
        {
            get => new SizeFloat(width, height);
            set { width = value.Width; height = value.Height; }
        }

        /// <summary>
        /// Gets the left coordinate.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public float Left => x;

        /// <summary>
        /// Gets the top coordinate.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public float Top => y;

        /// <summary>
        /// Gets the right coordinate (X + Width).
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public float Right => x + width;

        /// <summary>
        /// Gets the bottom coordinate (Y + Height).
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public float Bottom => y + height;

        /// <summary>
        /// Tests whether either dimension is zero or negative, regardless of location.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public bool IsEmpty => width <= 0 || height <= 0;
        #endregion

        #region Geometry
        /// <summary>
        /// Creates a rectangle from its left, top, right and bottom edges.
        /// </summary>
        /// <param name="left">Left coordinate.</param>
        /// <param name="top">Top coordinate.</param>
        /// <param name="right">Right coordinate.</param>
        /// <param name="bottom">Bottom coordinate.</param>
        /// <returns>Rectangle.</returns>
        public static RectangleFloat FromLTRB(float left, float top, float right, float bottom)
        {
            return new RectangleFloat(left, top, right - left, bottom - top);
        }

        /// <summary>
        /// Tests whether a point lies inside, excluding the right and bottom edges.
        /// </summary>
        /// <param name="x">Point X coordinate.</param>
        /// <param name="y">Point Y coordinate.</param>
        /// <returns>Boolean.</returns>
        public bool Contains(float x, float y)
        {
            return this.x <= x && x < Right && this.y <= y && y < Bottom;
        }

        /// <summary>
        /// Tests whether a point lies inside, excluding the right and bottom edges.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <returns>Boolean.</returns>
        public bool Contains(PointFloat pt) => Contains(pt.X, pt.Y);

        /// <summary>
        /// Tests whether the specified rectangle is entirely contained within this rectangle.
        /// </summary>
        /// <param name="rect">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public bool Contains(RectangleFloat rect)
        {
            return x <= rect.x && rect.Right <= Right
                && y <= rect.y && rect.Bottom <= Bottom;
        }

        /// <summary>
        /// Tests whether this rectangle overlaps the specified rectangle.
        /// Touching edges alone do not count as an overlap.
        /// </summary>
        /// <param name="rect">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public bool IntersectsWith(RectangleFloat rect)
        {
            return rect.x < Right && x < rect.Right
                && rect.y < Bottom && y < rect.Bottom;
        }

        /// <summary>
        /// Replaces this rectangle with its intersection with the specified rectangle.
        /// </summary>
        /// <param name="rect">Rectangle.</param>
        public void Intersect(RectangleFloat rect)
        {
            this = Intersect(rect, this);
        }

        /// <summary>
        /// Returns the intersection, or Empty if the rectangles are disjoint.
        /// Touching edges produce a zero-sized dimension at the touching location.
        /// </summary>
        /// <param name="a">First rectangle.</param>
        /// <param name="b">Second rectangle.</param>
        /// <returns>Intersection.</returns>
        public static RectangleFloat Intersect(RectangleFloat a, RectangleFloat b)
        {
            float left = Math.Max(a.Left, b.Left);
            float top = Math.Max(a.Top, b.Top);
            float right = Math.Min(a.Right, b.Right);
            float bottom = Math.Min(a.Bottom, b.Bottom);
            return right >= left && bottom >= top ? FromLTRB(left, top, right, bottom) : Empty;
        }

        /// <summary>
        /// Returns the smallest rectangle containing both specified rectangles.
        /// </summary>
        /// <param name="a">First rectangle.</param>
        /// <param name="b">Second rectangle.</param>
        /// <returns>Bounding rectangle.</returns>
        public static RectangleFloat Union(RectangleFloat a, RectangleFloat b)
        {
            return FromLTRB(Math.Min(a.Left, b.Left), Math.Min(a.Top, b.Top),
                Math.Max(a.Right, b.Right), Math.Max(a.Bottom, b.Bottom));
        }

        /// <summary>
        /// Expands this rectangle by the specified amount on each side.
        /// </summary>
        /// <param name="x">Horizontal amount per side.</param>
        /// <param name="y">Vertical amount per side.</param>
        public void Inflate(float x, float y)
        {
            this.x -= x;
            this.y -= y;
            width += 2 * x;
            height += 2 * y;
        }

        /// <summary>
        /// Expands this rectangle by the specified amount on each side.
        /// </summary>
        /// <param name="size">Horizontal and vertical amounts per side.</param>
        public void Inflate(SizeFloat size) => Inflate(size.Width, size.Height);

        /// <summary>
        /// Returns an expanded copy of the specified rectangle.
        /// </summary>
        /// <param name="rect">Rectangle.</param>
        /// <param name="x">Horizontal amount per side.</param>
        /// <param name="y">Vertical amount per side.</param>
        /// <returns>Expanded rectangle.</returns>
        public static RectangleFloat Inflate(RectangleFloat rect, float x, float y)
        {
            rect.Inflate(x, y);
            return rect;
        }

        /// <summary>
        /// Moves this rectangle by the specified displacement.
        /// </summary>
        /// <param name="x">Horizontal displacement.</param>
        /// <param name="y">Vertical displacement.</param>
        public void Offset(float x, float y)
        {
            this.x = this.x + x;
            this.y = this.y + y;
        }

        /// <summary>
        /// Moves this rectangle by the specified displacement.
        /// </summary>
        /// <param name="pos">Displacement.</param>
        public void Offset(PointFloat pos) => Offset(pos.X, pos.Y);
        #endregion

        #region Arithmetic
        /// <summary>
        /// Translates a rectangle by adding the specified offset to its position.
        /// </summary>
        /// <param name="point">Horizontal and vertical offsets to add.</param>
        /// <returns>Translated rectangle with the original width and height.</returns>
        public RectangleFloat Add(PointFloat point)
        {
            var rectangle = this;
            return new RectangleFloat
            {
                X = rectangle.X + point.X,
                Y = rectangle.Y + point.Y,
                Width = rectangle.Width,
                Height = rectangle.Height
            };
        }

        /// <summary>
        /// Translates a rectangle by subtracting the specified offset from its position.
        /// </summary>
        /// <param name="point">Horizontal and vertical offsets to subtract.</param>
        /// <returns>Translated rectangle with the original width and height.</returns>
        public RectangleFloat Sub(PointFloat point)
        {
            var rectangle = this;
            return new RectangleFloat
            {
                X = rectangle.X - point.X,
                Y = rectangle.Y - point.Y,
                Width = rectangle.Width,
                Height = rectangle.Height
            };
        }

        /// <summary>
        /// Translates each rectangle by adding the specified offset to its position.
        /// </summary>
        /// <param name="rectangles">Rectangles to translate.</param>
        /// <param name="point">Horizontal and vertical offsets to add.</param>
        /// <returns>New array of translated rectangles in the original order, with unchanged sizes.</returns>
        public static RectangleFloat[] Add(RectangleFloat[] rectangles, PointFloat point)
        {
            var count = rectangles.Length;
            var output = new RectangleFloat[count];

            for (int i = 0; i < count; i++)
            {
                output[i] = rectangles[i].Add(point);
            }

            return output;
        }

        /// <summary>
        /// Translates each rectangle by subtracting the specified offset from its position.
        /// </summary>
        /// <param name="rectangles">Rectangles to translate.</param>
        /// <param name="point">Horizontal and vertical offsets to subtract.</param>
        /// <returns>New array of translated rectangles in the original order, with unchanged sizes.</returns>
        public static RectangleFloat[] Sub(RectangleFloat[] rectangles, PointFloat point)
        {
            var count = rectangles.Length;
            var output = new RectangleFloat[count];

            for (int i = 0; i < count; i++)
            {
                output[i] = rectangles[i].Sub(point);
            }

            return output;
        }

        /// <summary>
        /// Returns the four corners of a rectangle.
        /// </summary>
        /// <returns>New array containing the top-left, top-right, bottom-right, and bottom-left corners, in that order.</returns>
        public PointFloat[] ToPoints()
        {
            var rectangle = this;
            return new PointFloat[]
            {
                new PointFloat (rectangle.Left, rectangle.Top),
                new PointFloat (rectangle.Right, rectangle.Top),
                new PointFloat (rectangle.Right, rectangle.Bottom),
                new PointFloat (rectangle.Left, rectangle.Bottom)
            };
        }

        /// <summary>
        /// Creates a rectangle from an array of four corner points.
        /// </summary>
        /// <param name="points">Four points ordered as top-left, top-right, bottom-right, and bottom-left.</param>
        /// <returns>Rectangle with its left and top edges taken from the first point and its right and bottom edges from the third point.</returns>
        /// <exception cref="ArgumentException">The array does not contain exactly four points.</exception>
        public static RectangleFloat FromPoints(PointFloat[] points)
        {
            if (points.Length != 4)
                throw new ArgumentException("A rectangle can only be built using four points");

            return RectangleFloat.FromLTRB(
                points[0].X,
                points[0].Y,
                points[2].X,
                points[2].Y);
        }

        /// <summary>
        /// Returns the location of a rectangle.
        /// </summary>
        /// <returns>Point containing the rectangle's X and Y coordinates.</returns>
        public PointFloat GetPoint()
        {
            var rectangle = this;
            return new PointFloat
            {
                X = rectangle.X,
                Y = rectangle.Y
            };
        }

        /// <summary>
        /// Calculates the product of a size's width and height.
        /// </summary>
        /// <param name="size">Size whose area is calculated.</param>
        /// <returns>Width multiplied by height, without taking the absolute value.</returns>
        public static float Area(SizeFloat size)
        {
            return size.Width * size.Height;
        }

        /// <summary>
        /// Calculates the product of a rectangle's width and height.
        /// </summary>
        /// <returns>Width multiplied by height, without taking the absolute value.</returns>
        public float Area()
        {
            var rectangle = this;
            return rectangle.Width * rectangle.Height;
        }

        /// <summary>
        /// Selects the rectangle with the largest width-height product.
        /// </summary>
        /// <param name="rectangles">Rectangles to compare by area.</param>
        /// <returns>First rectangle with the largest qualifying area, or the first input rectangle if none qualifies; <see cref="RectangleFloat.Empty"/> for an empty array.</returns>
        public static RectangleFloat Max(params RectangleFloat[] rectangles)
        {
            var length = rectangles.Length;
            var rectangle = RectangleFloat.Empty;
            var area = float.MinValue;
            var max = 0;

            for (int i = 0; i < length; i++)
            {
                rectangle = rectangles[i];

                if (rectangle.IsEmpty)
                    continue;

                var current = rectangle.Area();

                if (current > area)
                {
                    max = i;
                    area = current;
                }
            }

            return length > 0 ? rectangles[max] : rectangle;
        }

        /// <summary>
        /// Selects the rectangle with the smallest width-height product.
        /// </summary>
        /// <param name="rectangles">Rectangles to compare by area.</param>
        /// <returns>First rectangle with the smallest qualifying area, or <see cref="RectangleFloat.Empty"/> if none qualifies.</returns>
        public static RectangleFloat Min(params RectangleFloat[] rectangles)
        {
            var length = rectangles.Length;
            var rectangle = RectangleFloat.Empty;
            var area = float.MaxValue;
            var min = -1;

            for (int i = 0; i < length; i++)
            {
                rectangle = rectangles[i];

                if (rectangle.IsEmpty)
                    continue;

                var current = rectangle.Area();

                if (current < area)
                {
                    min = i;
                    area = current;
                }
            }

            return min >= 0 ? rectangles[min] : RectangleInt.Empty;
        }

        /// <summary>
        /// Converts a rectangle to a square using its larger dimension.
        /// </summary>
        /// <returns>Square with a side equal to the larger of the original width and height.</returns>
        public RectangleFloat ToBox()
        {
            var rectangle = this;
            var max = Math.Max(rectangle.Width, rectangle.Height);
            var dx = max - rectangle.Width;
            var dy = max - rectangle.Height;

            return new RectangleFloat
            {
                X = rectangle.X - dx / 2,
                Y = rectangle.Y - dy / 2,
                Width = rectangle.Width + dx,
                Height = rectangle.Height + dy
            };
        }

        /// <summary>
        /// Resizes a rectangle by a relative change in both dimensions.
        /// </summary>
        /// <param name="scale">Relative change in width and height; for example, 0.1 increases both by 10% before truncation.</param>
        /// <returns>Resized rectangle with coordinates and dimensions truncated toward zero to integers.</returns>
        public RectangleFloat ToBox(float scale)
        {
            var rectangle = this;
            float gainX = rectangle.Width * scale;
            float gainY = rectangle.Height * scale;

            return new RectangleFloat(
                (int)(rectangle.X - gainX / 2),
                (int)(rectangle.Y - gainY / 2),
                (int)(rectangle.Width + gainX),
                (int)(rectangle.Height + gainY)
                );
        }

        /// <summary>
        /// Converts each rectangle to a square using its larger dimension.
        /// </summary>
        /// <param name="rectangles">Rectangles to convert.</param>
        /// <returns>New array of squares in the original order, each with a side equal to the larger original dimension.</returns>
        public static RectangleFloat[] ToBox(params RectangleFloat[] rectangles)
        {
            int length = rectangles.Length;
            var newRectangles = new RectangleFloat[length];

            for (int i = 0; i < length; i++)
            {
                newRectangles[i] = rectangles[i].ToBox();
            }

            return newRectangles;
        }

        /// <summary>
        /// Resizes each rectangle by a relative change in both dimensions.
        /// </summary>
        /// <param name="factor">Relative change in width and height; for example, 0.1 increases both by 10% before truncation.</param>
        /// <param name="rectangles">Rectangles to resize.</param>
        /// <returns>New array of resized rectangles in the original order.</returns>
        public static RectangleFloat[] ToBox(float factor, params RectangleFloat[] rectangles)
        {
            int length = rectangles.Length;
            var newRectangles = new RectangleFloat[length];

            for (int i = 0; i < length; i++)
            {
                newRectangles[i] = rectangles[i].ToBox(factor);
            }

            return newRectangles;
        }

        /// <summary>
        /// Calculates the intersection-over-union (IoU) ratio of two rectangles.
        /// </summary>
        /// <param name="b">Second rectangle, with nonnegative width and height.</param>
        /// <returns>Intersection area divided by union area, or zero if the intersection area is zero.</returns>
        public float IoU(RectangleFloat b)
        {
            var a = this;
            var xA = Math.Max(a.Left, b.Left);
            var yA = Math.Max(a.Top, b.Top);
            var xB = Math.Min(a.Right, b.Right);
            var yB = Math.Min(a.Bottom, b.Bottom);

            var interArea = Math.Abs(Math.Max(xB - xA, 0) * (float)Math.Max(yB - yA, 0));

            if (interArea == 0)
                return 0;

            var boxAArea = Math.Abs((a.Right - a.Left) * (float)(a.Bottom - a.Top));
            var boxBArea = Math.Abs((b.Right - b.Left) * (float)(b.Bottom - b.Top));

            return interArea / (float)(boxAArea + boxBArea - interArea);
        }

        /// <summary>
        /// Resizes a rectangle by independent relative changes in width and height.
        /// </summary>
        /// <param name="kx">Relative width change; for example, 0.1 adds 10% of the width before truncation. Defaults to zero.</param>
        /// <param name="ky">Relative height change; for example, 0.1 adds 10% of the height before truncation. Defaults to zero.</param>
        /// <returns>Rectangle with the computed dimension changes added and half of each change subtracted from its position.</returns>
        public RectangleFloat Scale(float kx = 0.0f, float ky = 0.0f)
        {
            var rectangle = this;
            var x = rectangle.X;
            var y = rectangle.Y;
            var w = rectangle.Width;
            var h = rectangle.Height;

            var dw = (int)(w * kx);
            var dh = (int)(h * ky);

            return new RectangleFloat
            {
                X = x - dw / 2,
                Y = y - dh / 2,
                Width = w + dw,
                Height = h + dh,
            };
        }

        /// <summary>
        /// Converts a rectangle to a square using its diagonal length.
        /// </summary>
        /// <returns>Square with a side equal to the original diagonal length truncated toward zero to an integer.</returns>
        public RectangleFloat Scale()
        {
            var rectangle = this;
            var r = (int)Math.Sqrt(rectangle.Width * rectangle.Width + rectangle.Height * rectangle.Height);
            var dx = r - rectangle.Width;
            var dy = r - rectangle.Height;

            var x = rectangle.X - dx / 2;
            var y = rectangle.Y - dy / 2;
            var w = rectangle.Width + dx;
            var h = rectangle.Height + dy;

            return new RectangleFloat
            {
                X = x,
                Y = y,
                Width = w,
                Height = h
            };
        }

        /// <summary>
        /// Clips a rectangle against the bounds of another rectangle.
        /// </summary>
        /// <param name="second">Clipping bounds, expected to have nonnegative width and height.</param>
        /// <returns>Rectangle starting at the maximum left and top coordinates, with each intersection dimension clamped to zero.</returns>
        public RectangleFloat Clamp(RectangleFloat second)
        {
            var first = this;
            if (first.Width < 0) { first.X += first.Width; first.Width = -first.Width; }
            if (first.Height < 0) { first.Y += first.Height; first.Height = -first.Height; }

            float x = Math.Max(first.X, second.Left);
            float y = Math.Max(first.Y, second.Top);

            float right = Math.Min(first.Right, second.Right);
            float bottom = Math.Min(first.Bottom, second.Bottom);

            float w = Math.Max(0, right - x);
            float h = Math.Max(0, bottom - y);

            return new RectangleFloat(x, y, w, h);
        }
        #endregion

        #region Conversions
        /// <summary>
        /// Initializes a rectangle from vector components X, Y, Z (width), W (height).
        /// </summary>
        /// <param name="vector">Location and dimensions.</param>
        public RectangleFloat(Vector4 vector) : this(vector.X, vector.Y, vector.Z, vector.W)
        {
        }

        /// <summary>
        /// Returns a vector containing X, Y, Width and Height.
        /// </summary>
        /// <returns>Four-component vector.</returns>
        public Vector4 ToVector4() => new Vector4(x, y, width, height);

        /// <summary>
        /// Converts a rectangle to a vector.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        public static explicit operator Vector4(RectangleFloat rectangle) => rectangle.ToVector4();

        /// <summary>
        /// Converts a vector to a rectangle.
        /// </summary>
        /// <param name="vector">Location and dimensions.</param>
        public static explicit operator RectangleFloat(Vector4 vector) => new RectangleFloat(vector);

        /// <summary>
        /// Converts integer components to single precision.
        /// </summary>
        /// <param name="r">Integer rectangle.</param>
        public static implicit operator RectangleFloat(RectangleInt r)
        {
            return new RectangleFloat(r.X, r.Y, r.Width, r.Height);
        }

        #endregion

        #region Equality and overrides
        /// <summary>
        /// Compares all four components.
        /// </summary>
        /// <param name="left">First rectangle.</param>
        /// <param name="right">Second rectangle.</param>
        /// <returns>Boolean.</returns>
        public static bool operator ==(RectangleFloat left, RectangleFloat right)
        {
            return left.x == right.x && left.y == right.y && left.width == right.width && left.height == right.height;
        }

        /// <summary>
        /// Tests whether any component differs.
        /// </summary>
        /// <param name="left">First rectangle.</param>
        /// <param name="right">Second rectangle.</param>
        /// <returns>Boolean.</returns>
        public static bool operator !=(RectangleFloat left, RectangleFloat right) => !(left == right);

        /// <summary>
        /// Tests whether another rectangle has the same components.
        /// </summary>
        /// <param name="other">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public bool Equals(RectangleFloat other) => this == other;

        /// <summary>
        /// Tests whether an object is an equal RectangleFloat.
        /// </summary>
        /// <param name="obj">Object.</param>
        /// <returns>Boolean.</returns>
        public override bool Equals(object obj) => obj is RectangleFloat other && Equals(other);

        /// <summary>
        /// Returns a hash code based on all four components.
        /// </summary>
        /// <returns>Hash code.</returns>
        public override int GetHashCode()
        {
            unchecked
            {
                int hash = x.GetHashCode();
                hash = (hash * 397) ^ y.GetHashCode();
                hash = (hash * 397) ^ width.GetHashCode();
                return (hash * 397) ^ height.GetHashCode();
            }
        }

        /// <summary>
        /// Returns a string containing the location and dimensions.
        /// </summary>
        /// <returns>Text representation.</returns>
        public override string ToString()
        {
            return string.Format("{{X={0},Y={1},Width={2},Height={3}}}", x, y, width, height);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of RectangleFloat.
        /// </summary>
        /// <returns>Rectangle copy.</returns>
        object ICloneable.Clone()
        {
            return new RectangleFloat(x, y, width, height);
        }

        /// <summary>
        /// Creates a copy of RectangleFloat.
        /// </summary>
        /// <returns>Rectangle copy.</returns>
        public RectangleFloat Clone()
        {
            return new RectangleFloat(x, y, width, height);
        }

        #endregion
    }
}
