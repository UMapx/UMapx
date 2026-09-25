using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a rectangle with integer coordinates and dimensions.
    /// </summary>
    [Serializable]
    public struct RectangleInt : IEquatable<RectangleInt>, ICloneable
    {
        #region Private data
        private int x;
        private int y;
        private int width;
        private int height;
        #endregion

        #region Structure components
        /// <summary>
        /// Represents a rectangle whose coordinates and dimensions are zero.
        /// </summary>
        public static readonly RectangleInt Empty;

        /// <summary>
        /// Initializes a rectangle with the specified location and size.
        /// </summary>
        /// <param name="x">Left coordinate.</param>
        /// <param name="y">Top coordinate.</param>
        /// <param name="width">Width.</param>
        /// <param name="height">Height.</param>
        public RectangleInt(int x, int y, int width, int height)
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
        public RectangleInt(PointInt location, SizeInt size)
            : this(location.X, location.Y, size.Width, size.Height)
        {
        }

        /// <summary>
        /// Gets or sets the left coordinate.
        /// </summary>
        public int X { get => x; set => x = value; }

        /// <summary>
        /// Gets or sets the top coordinate.
        /// </summary>
        public int Y { get => y; set => y = value; }

        /// <summary>
        /// Gets or sets the width.
        /// </summary>
        public int Width { get => width; set => width = value; }

        /// <summary>
        /// Gets or sets the height.
        /// </summary>
        public int Height { get => height; set => height = value; }

        /// <summary>
        /// Gets or sets the upper-left corner.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public PointInt Location
        {
            get => new PointInt(x, y);
            set { x = value.X; y = value.Y; }
        }

        /// <summary>
        /// Gets or sets the dimensions.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public SizeInt Size
        {
            get => new SizeInt(width, height);
            set { width = value.Width; height = value.Height; }
        }

        /// <summary>
        /// Gets the left coordinate.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public int Left => x;

        /// <summary>
        /// Gets the top coordinate.
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public int Top => y;

        /// <summary>
        /// Gets the right coordinate (X + Width).
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public int Right => unchecked(x + width);

        /// <summary>
        /// Gets the bottom coordinate (Y + Height).
        /// </summary>
        [System.ComponentModel.Browsable(false)]
        public int Bottom => unchecked(y + height);

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
        public static RectangleInt FromLTRB(int left, int top, int right, int bottom)
        {
            return new RectangleInt(left, top, unchecked(right - left), unchecked(bottom - top));
        }

        /// <summary>
        /// Tests whether a point lies inside, excluding the right and bottom edges.
        /// </summary>
        /// <param name="x">Point X coordinate.</param>
        /// <param name="y">Point Y coordinate.</param>
        /// <returns>Boolean.</returns>
        public bool Contains(int x, int y)
        {
            return this.x <= x && x < Right && this.y <= y && y < Bottom;
        }

        /// <summary>
        /// Tests whether a point lies inside, excluding the right and bottom edges.
        /// </summary>
        /// <param name="pt">Point.</param>
        /// <returns>Boolean.</returns>
        public bool Contains(PointInt pt) => Contains(pt.X, pt.Y);

        /// <summary>
        /// Tests whether the specified rectangle is entirely contained within this rectangle.
        /// </summary>
        /// <param name="rect">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public bool Contains(RectangleInt rect)
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
        public bool IntersectsWith(RectangleInt rect)
        {
            return rect.x < Right && x < rect.Right
                && rect.y < Bottom && y < rect.Bottom;
        }

        /// <summary>
        /// Replaces this rectangle with its intersection with the specified rectangle.
        /// </summary>
        /// <param name="rect">Rectangle.</param>
        public void Intersect(RectangleInt rect)
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
        public static RectangleInt Intersect(RectangleInt a, RectangleInt b)
        {
            int left = Math.Max(a.Left, b.Left);
            int top = Math.Max(a.Top, b.Top);
            int right = Math.Min(a.Right, b.Right);
            int bottom = Math.Min(a.Bottom, b.Bottom);
            return right >= left && bottom >= top ? FromLTRB(left, top, right, bottom) : Empty;
        }

        /// <summary>
        /// Returns the smallest rectangle containing both specified rectangles.
        /// </summary>
        /// <param name="a">First rectangle.</param>
        /// <param name="b">Second rectangle.</param>
        /// <returns>Bounding rectangle.</returns>
        public static RectangleInt Union(RectangleInt a, RectangleInt b)
        {
            return FromLTRB(Math.Min(a.Left, b.Left), Math.Min(a.Top, b.Top),
                Math.Max(a.Right, b.Right), Math.Max(a.Bottom, b.Bottom));
        }

        /// <summary>
        /// Expands this rectangle by the specified amount on each side.
        /// </summary>
        /// <param name="width">Horizontal amount per side.</param>
        /// <param name="height">Vertical amount per side.</param>
        public void Inflate(int width, int height)
        {
            unchecked
            {
                x -= width;
                y -= height;
                this.width += 2 * width;
                this.height += 2 * height;
            }
        }

        /// <summary>
        /// Expands this rectangle by the specified amount on each side.
        /// </summary>
        /// <param name="size">Horizontal and vertical amounts per side.</param>
        public void Inflate(SizeInt size) => Inflate(size.Width, size.Height);

        /// <summary>
        /// Returns an expanded copy of the specified rectangle.
        /// </summary>
        /// <param name="rect">Rectangle.</param>
        /// <param name="x">Horizontal amount per side.</param>
        /// <param name="y">Vertical amount per side.</param>
        /// <returns>Expanded rectangle.</returns>
        public static RectangleInt Inflate(RectangleInt rect, int x, int y)
        {
            rect.Inflate(x, y);
            return rect;
        }

        /// <summary>
        /// Moves this rectangle by the specified displacement.
        /// </summary>
        /// <param name="x">Horizontal displacement.</param>
        /// <param name="y">Vertical displacement.</param>
        public void Offset(int x, int y)
        {
            this.x = unchecked(this.x + x);
            this.y = unchecked(this.y + y);
        }

        /// <summary>
        /// Moves this rectangle by the specified displacement.
        /// </summary>
        /// <param name="pos">Displacement.</param>
        public void Offset(PointInt pos) => Offset(pos.X, pos.Y);
        #endregion

        #region Arithmetic
        /// <summary>
        /// Translates a rectangle by adding the specified offset to its position.
        /// </summary>
        /// <param name="point">Horizontal and vertical offsets to add.</param>
        /// <returns>Translated rectangle with the original width and height.</returns>
        public RectangleInt Add(PointInt point)
        {
            var rectangle = this;
            return new RectangleInt
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
        public RectangleInt Sub(PointInt point)
        {
            var rectangle = this;
            return new RectangleInt
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
        public static RectangleInt[] Add(RectangleInt[] rectangles, PointInt point)
        {
            var count = rectangles.Length;
            var output = new RectangleInt[count];

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
        public static RectangleInt[] Sub(RectangleInt[] rectangles, PointInt point)
        {
            var count = rectangles.Length;
            var output = new RectangleInt[count];

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
        public PointInt[] ToPoints()
        {
            var rectangle = this;
            return new PointInt[]
            {
                new PointInt (rectangle.Left, rectangle.Top),
                new PointInt (rectangle.Right, rectangle.Top),
                new PointInt (rectangle.Right, rectangle.Bottom),
                new PointInt (rectangle.Left, rectangle.Bottom)
            };
        }

        /// <summary>
        /// Creates a rectangle from an array of four corner points.
        /// </summary>
        /// <param name="points">Four points ordered as top-left, top-right, bottom-right, and bottom-left.</param>
        /// <returns>RectangleInt with its left and top edges taken from the first point and its right and bottom edges from the third point.</returns>
        /// <exception cref="ArgumentException">The array does not contain exactly four points.</exception>
        public static RectangleInt FromPoints(PointInt[] points)
        {
            if (points.Length != 4)
                throw new ArgumentException("A rectangle can only be built using four points");

            return RectangleInt.FromLTRB(
                points[0].X,
                points[0].Y,
                points[2].X,
                points[2].Y);
        }

        /// <summary>
        /// Returns the location of a rectangle.
        /// </summary>
        /// <returns>PointInt containing the rectangle's X and Y coordinates.</returns>
        public PointInt GetPoint()
        {
            var rectangle = this;
            return new PointInt
            {
                X = rectangle.X,
                Y = rectangle.Y
            };
        }

        /// <summary>
        /// Calculates the product of a size's width and height.
        /// </summary>
        /// <param name="size">SizeInt whose area is calculated.</param>
        /// <returns>Width multiplied by height, without taking the absolute value.</returns>
        public static int Area(SizeInt size)
        {
            return size.Width * size.Height;
        }

        /// <summary>
        /// Calculates the product of a rectangle's width and height.
        /// </summary>
        /// <returns>Width multiplied by height, without taking the absolute value.</returns>
        public int Area()
        {
            var rectangle = this;
            return rectangle.Width * rectangle.Height;
        }

        /// <summary>
        /// Selects the rectangle with the largest width-height product.
        /// </summary>
        /// <param name="rectangles">Rectangles to compare by area.</param>
        /// <returns>First rectangle with the largest qualifying area, or the first input rectangle if none qualifies; <see cref="RectangleInt.Empty"/> for an empty array.</returns>
        public static RectangleInt Max(params RectangleInt[] rectangles)
        {
            var length = rectangles.Length;
            var rectangle = RectangleInt.Empty;
            var area = int.MinValue;
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
        /// <returns>First rectangle with the smallest qualifying area, or <see cref="RectangleInt.Empty"/> if none qualifies.</returns>
        public static RectangleInt Min(params RectangleInt[] rectangles)
        {
            var length = rectangles.Length;
            var rectangle = RectangleInt.Empty;
            var area = int.MaxValue;
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
        public RectangleInt ToBox()
        {
            var rectangle = this;
            var max = Math.Max(rectangle.Width, rectangle.Height);
            var dx = max - rectangle.Width;
            var dy = max - rectangle.Height;

            return new RectangleInt
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
        public RectangleInt ToBox(float scale)
        {
            var rectangle = this;
            float gainX = rectangle.Width * scale;
            float gainY = rectangle.Height * scale;

            return new RectangleInt(
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
        public static RectangleInt[] ToBox(params RectangleInt[] rectangles)
        {
            int length = rectangles.Length;
            var newRectangles = new RectangleInt[length];

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
        public static RectangleInt[] ToBox(float factor, params RectangleInt[] rectangles)
        {
            int length = rectangles.Length;
            var newRectangles = new RectangleInt[length];

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
        public float IoU(RectangleInt b)
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
        /// <returns>RectangleInt with the computed dimension changes added and half of each change subtracted from its position.</returns>
        public RectangleInt Scale(float kx = 0.0f, float ky = 0.0f)
        {
            var rectangle = this;
            var x = rectangle.X;
            var y = rectangle.Y;
            var w = rectangle.Width;
            var h = rectangle.Height;

            var dw = (int)(w * kx);
            var dh = (int)(h * ky);

            return new RectangleInt
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
        public RectangleInt Scale()
        {
            var rectangle = this;
            var r = (int)Math.Sqrt(rectangle.Width * rectangle.Width + rectangle.Height * rectangle.Height);
            var dx = r - rectangle.Width;
            var dy = r - rectangle.Height;

            var x = rectangle.X - dx / 2;
            var y = rectangle.Y - dy / 2;
            var w = rectangle.Width + dx;
            var h = rectangle.Height + dy;

            return new RectangleInt
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
        /// <returns>RectangleInt starting at the maximum left and top coordinates, with each intersection dimension clamped to zero.</returns>
        public RectangleInt Clamp(RectangleInt second)
        {
            var first = this;
            if (first.Width < 0) { first.X += first.Width; first.Width = -first.Width; }
            if (first.Height < 0) { first.Y += first.Height; first.Height = -first.Height; }

            int x = Math.Max(first.X, second.Left);
            int y = Math.Max(first.Y, second.Top);

            int right = Math.Min(first.Right, second.Right);
            int bottom = Math.Min(first.Bottom, second.Bottom);

            int w = Math.Max(0, right - x);
            int h = Math.Max(0, bottom - y);

            return new RectangleInt(x, y, w, h);
        }
        #endregion

        #region Conversions
        /// <summary>
        /// Rounds each floating-point component toward positive infinity.
        /// </summary>
        /// <param name="value">Floating-point rectangle.</param>
        /// <returns>Integer rectangle.</returns>
        public static RectangleInt Ceiling(RectangleFloat value)
        {
            return new RectangleInt(unchecked((int)Math.Ceiling(value.X)),
                unchecked((int)Math.Ceiling(value.Y)),
                unchecked((int)Math.Ceiling(value.Width)),
                unchecked((int)Math.Ceiling(value.Height)));
        }

        /// <summary>
        /// Rounds each floating-point component to the nearest integer, with midpoint ties to even.
        /// </summary>
        /// <param name="value">Floating-point rectangle.</param>
        /// <returns>Integer rectangle.</returns>
        public static RectangleInt Round(RectangleFloat value)
        {
            return new RectangleInt(unchecked((int)Math.Round(value.X)),
                unchecked((int)Math.Round(value.Y)),
                unchecked((int)Math.Round(value.Width)),
                unchecked((int)Math.Round(value.Height)));
        }

        /// <summary>
        /// Truncates each floating-point component toward zero.
        /// </summary>
        /// <param name="value">Floating-point rectangle.</param>
        /// <returns>Integer rectangle.</returns>
        public static RectangleInt Truncate(RectangleFloat value)
        {
            return new RectangleInt(unchecked((int)value.X), unchecked((int)value.Y),
                unchecked((int)value.Width), unchecked((int)value.Height));
        }

        #endregion

        #region Equality and overrides
        /// <summary>
        /// Compares all four components.
        /// </summary>
        /// <param name="left">First rectangle.</param>
        /// <param name="right">Second rectangle.</param>
        /// <returns>Boolean.</returns>
        public static bool operator ==(RectangleInt left, RectangleInt right)
        {
            return left.x == right.x && left.y == right.y && left.width == right.width && left.height == right.height;
        }

        /// <summary>
        /// Tests whether any component differs.
        /// </summary>
        /// <param name="left">First rectangle.</param>
        /// <param name="right">Second rectangle.</param>
        /// <returns>Boolean.</returns>
        public static bool operator !=(RectangleInt left, RectangleInt right) => !(left == right);

        /// <summary>
        /// Tests whether another rectangle has the same components.
        /// </summary>
        /// <param name="other">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public bool Equals(RectangleInt other) => this == other;

        /// <summary>
        /// Tests whether an object is an equal RectangleInt.
        /// </summary>
        /// <param name="obj">Object.</param>
        /// <returns>Boolean.</returns>
        public override bool Equals(object obj) => obj is RectangleInt other && Equals(other);

        /// <summary>
        /// Returns a hash code based on all four components.
        /// </summary>
        /// <returns>Hash code.</returns>
        public override int GetHashCode()
        {
            return new System.Drawing.Rectangle(x, y, width, height).GetHashCode();
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
        /// Creates a copy of RectangleInt.
        /// </summary>
        /// <returns>Rectangle copy.</returns>
        object ICloneable.Clone()
        {
            return new RectangleInt(x, y, width, height);
        }

        /// <summary>
        /// Creates a copy of RectangleInt.
        /// </summary>
        /// <returns>Rectangle copy.</returns>
        public RectangleInt Clone()
        {
            return new RectangleInt(x, y, width, height);
        }

        #endregion
    }
}
