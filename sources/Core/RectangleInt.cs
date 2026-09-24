using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a rectangle with integer coordinates and dimensions.
    /// </summary>
    /// <remarks>
    /// Standard geometry members follow System.Drawing.Rectangle semantics, except IsEmpty:
    /// both rectangle types treat zero or negative dimensions as empty.
    /// Additional arithmetic members return new values without changing this rectangle.
    /// </remarks>
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
        public int X { readonly get => x; set => x = value; }

        /// <summary>
        /// Gets or sets the top coordinate.
        /// </summary>
        public int Y { readonly get => y; set => y = value; }

        /// <summary>
        /// Gets or sets the width.
        /// </summary>
        public int Width { readonly get => width; set => width = value; }

        /// <summary>
        /// Gets or sets the height.
        /// </summary>
        public int Height { readonly get => height; set => height = value; }

        /// <summary>
        /// Gets or sets the upper-left corner.
        /// </summary>
        public PointInt Location
        {
            readonly get => new PointInt(x, y);
            set { x = value.X; y = value.Y; }
        }

        /// <summary>
        /// Gets or sets the dimensions.
        /// </summary>
        public SizeInt Size
        {
            readonly get => new SizeInt(width, height);
            set { width = value.Width; height = value.Height; }
        }

        /// <summary>
        /// Gets the left coordinate.
        /// </summary>
        public readonly int Left => x;

        /// <summary>
        /// Gets the top coordinate.
        /// </summary>
        public readonly int Top => y;

        /// <summary>
        /// Gets the right coordinate (X + Width).
        /// </summary>
        public readonly int Right => unchecked(x + width);

        /// <summary>
        /// Gets the bottom coordinate (Y + Height).
        /// </summary>
        public readonly int Bottom => unchecked(y + height);

        /// <summary>
        /// Tests whether either dimension is zero or negative, regardless of location.
        /// </summary>
        /// <remarks>
        /// Uses the same rule as RectangleFloat. Compare with Empty to test whether all components are zero.
        /// </remarks>
        public readonly bool IsEmpty => width <= 0 || height <= 0;
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
        public readonly bool Contains(int x, int y)
        {
            return this.x <= x && x < Right && this.y <= y && y < Bottom;
        }

        /// <summary>
        /// Tests whether a point lies inside, excluding the right and bottom edges.
        /// </summary>
        /// <param name="point">Point.</param>
        /// <returns>Boolean.</returns>
        public readonly bool Contains(PointInt point) => Contains(point.X, point.Y);

        /// <summary>
        /// Tests whether the specified rectangle is entirely contained within this rectangle.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public readonly bool Contains(RectangleInt rectangle)
        {
            return x <= rectangle.x && rectangle.Right <= Right
                && y <= rectangle.y && rectangle.Bottom <= Bottom;
        }

        /// <summary>
        /// Tests whether this rectangle overlaps the specified rectangle.
        /// Touching edges alone do not count as an overlap.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public readonly bool IntersectsWith(RectangleInt rectangle)
        {
            return rectangle.x < Right && x < rectangle.Right
                && rectangle.y < Bottom && y < rectangle.Bottom;
        }

        /// <summary>
        /// Replaces this rectangle with its intersection with the specified rectangle.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        public void Intersect(RectangleInt rectangle)
        {
            this = Intersect(rectangle, this);
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
        /// <param name="x">Horizontal amount per side.</param>
        /// <param name="y">Vertical amount per side.</param>
        public void Inflate(int x, int y)
        {
            unchecked
            {
                this.x -= x;
                this.y -= y;
                width += 2 * x;
                height += 2 * y;
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
        /// <param name="rectangle">Rectangle.</param>
        /// <param name="x">Horizontal amount per side.</param>
        /// <param name="y">Vertical amount per side.</param>
        /// <returns>Expanded rectangle.</returns>
        public static RectangleInt Inflate(RectangleInt rectangle, int x, int y)
        {
            rectangle.Inflate(x, y);
            return rectangle;
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
        /// <param name="point">Displacement.</param>
        public void Offset(PointInt point) => Offset(point.X, point.Y);
        #endregion

        #region Arithmetic
        /// <summary>
        /// Returns a translated copy of this rectangle.
        /// </summary>
        /// <param name="point">Displacement to add.</param>
        /// <returns>Translated rectangle.</returns>
        public readonly RectangleInt Add(PointInt point)
        {
            return new RectangleInt(unchecked(x + point.X), unchecked(y + point.Y), width, height);
        }

        /// <summary>
        /// Returns a copy translated by the negative of the displacement.
        /// </summary>
        /// <param name="point">Displacement to subtract.</param>
        /// <returns>Translated rectangle.</returns>
        public readonly RectangleInt Sub(PointInt point)
        {
            return new RectangleInt(unchecked(x - point.X), unchecked(y - point.Y), width, height);
        }

        /// <summary>
        /// Adds a displacement without changing the original rectangle.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        /// <param name="point">Displacement.</param>
        /// <returns>Translated rectangle.</returns>
        public static RectangleInt operator +(RectangleInt rectangle, PointInt point) => rectangle.Add(point);

        /// <summary>
        /// Subtracts a displacement without changing the original rectangle.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        /// <param name="point">Displacement.</param>
        /// <returns>Translated rectangle.</returns>
        public static RectangleInt operator -(RectangleInt rectangle, PointInt point) => rectangle.Sub(point);

        /// <summary>
        /// Returns translated copies of an array of rectangles.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <param name="point">Displacement to add.</param>
        /// <returns>New array of translated rectangles.</returns>
        public static RectangleInt[] Add(RectangleInt[] rectangles, PointInt point)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = new RectangleInt[rectangles.Length];
            for (int i = 0; i < result.Length; i++) result[i] = rectangles[i].Add(point);
            return result;
        }

        /// <summary>
        /// Returns copies of an array translated by the negative displacement.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <param name="point">Displacement to subtract.</param>
        /// <returns>New array of translated rectangles.</returns>
        public static RectangleInt[] Sub(RectangleInt[] rectangles, PointInt point)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = new RectangleInt[rectangles.Length];
            for (int i = 0; i < result.Length; i++) result[i] = rectangles[i].Sub(point);
            return result;
        }

        /// <summary>
        /// Returns corners clockwise: top-left, top-right, bottom-right, bottom-left.
        /// </summary>
        /// <returns>Four corner points.</returns>
        public readonly PointInt[] ToPoints()
        {
            return new[] { new PointInt(Left, Top), new PointInt(Right, Top),
                new PointInt(Right, Bottom), new PointInt(Left, Bottom) };
        }

        /// <summary>
        /// Creates a rectangle from four ordered corners, using corners zero and two.
        /// </summary>
        /// <param name="points">Four corners in the order returned by ToPoints.</param>
        /// <returns>Rectangle.</returns>
        /// <exception cref="ArgumentNullException">The points array is null.</exception>
        /// <exception cref="ArgumentException">The array does not contain four points.</exception>
        public static RectangleInt FromPoints(PointInt[] points)
        {
            if (points == null) throw new ArgumentNullException(nameof(points));
            if (points.Length != 4)
                throw new ArgumentException("A rectangle can only be built using four points.", nameof(points));
            return FromLTRB(points[0].X, points[0].Y, points[2].X, points[2].Y);
        }

        /// <summary>
        /// Returns the upper-left corner.
        /// </summary>
        /// <returns>Point.</returns>
        public readonly PointInt GetPoint() => Location;

        /// <summary>
        /// Returns Width multiplied by Height, using unchecked integer arithmetic.
        /// </summary>
        /// <returns>Signed area.</returns>
        public readonly int Area() => unchecked(width * height);

        /// <summary>
        /// Returns the product of the specified dimensions.
        /// </summary>
        /// <param name="size">Dimensions.</param>
        /// <returns>Signed area, using unchecked integer arithmetic.</returns>
        public static int Area(SizeInt size) => unchecked(size.Width * size.Height);

        /// <summary>
        /// Returns the first non-empty rectangle with the largest signed area.
        /// Products are compared in 64-bit arithmetic to avoid area overflow.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <returns>Selected rectangle, or Empty if no non-empty rectangle exists.</returns>
        public static RectangleInt Max(params RectangleInt[] rectangles)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = Empty;
            long bestArea = 0;
            bool found = false;
            foreach (var rectangle in rectangles)
            {
                if (rectangle.IsEmpty) continue;
                long area = (long)rectangle.width * rectangle.height;
                if (!found || area > bestArea)
                {
                    result = rectangle;
                    bestArea = area;
                    found = true;
                }
            }
            return result;
        }

        /// <summary>
        /// Returns the first non-empty rectangle with the smallest signed area.
        /// Products are compared in 64-bit arithmetic to avoid area overflow.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <returns>Selected rectangle, or Empty if no non-empty rectangle exists.</returns>
        public static RectangleInt Min(params RectangleInt[] rectangles)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = Empty;
            long bestArea = 0;
            bool found = false;
            foreach (var rectangle in rectangles)
            {
                if (rectangle.IsEmpty) continue;
                long area = (long)rectangle.width * rectangle.height;
                if (!found || area < bestArea)
                {
                    result = rectangle;
                    bestArea = area;
                    found = true;
                }
            }
            return result;
        }

        /// <summary>
        /// Expands the shorter dimension to form a square around the same center.
        /// Integer half-displacements are truncated toward zero.
        /// </summary>
        /// <returns>Square rectangle.</returns>
        public readonly RectangleInt ToBox()
        {
            int side = Math.Max(width, height);
            return new RectangleInt(x - (side - width) / 2, y - (side - height) / 2, side, side);
        }

        /// <summary>
        /// Grows each dimension by the specified fraction around the center.
        /// A scale of zero preserves the rectangle; one doubles its dimensions.
        /// Each resulting component is truncated toward zero; this overload does not form a square.
        /// </summary>
        /// <param name="scale">Fractional increase in both dimensions.</param>
        /// <returns>Resized rectangle.</returns>
        public readonly RectangleInt ToBox(float scale)
        {
            float dx = width * scale;
            float dy = height * scale;
            return new RectangleInt((int)(x - dx / 2), (int)(y - dy / 2), (int)(width + dx), (int)(height + dy));
        }

        /// <summary>
        /// Returns square copies of the specified rectangles.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <returns>New array of square rectangles.</returns>
        public static RectangleInt[] ToBox(params RectangleInt[] rectangles)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = new RectangleInt[rectangles.Length];
            for (int i = 0; i < result.Length; i++) result[i] = rectangles[i].ToBox();
            return result;
        }

        /// <summary>
        /// Returns copies grown by the specified fraction around their centers.
        /// </summary>
        /// <param name="factor">Fractional increase in both dimensions.</param>
        /// <param name="rectangles">Rectangles.</param>
        /// <returns>New array of resized rectangles.</returns>
        public static RectangleInt[] ToBox(float factor, params RectangleInt[] rectangles)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = new RectangleInt[rectangles.Length];
            for (int i = 0; i < result.Length; i++) result[i] = rectangles[i].ToBox(factor);
            return result;
        }

        /// <summary>
        /// Returns intersection area divided by union area.
        /// Degenerate or disjoint rectangles return zero. Intermediate calculations
        /// use double precision to avoid coordinate and area overflow.
        /// </summary>
        /// <param name="rectangle">Other rectangle.</param>
        /// <returns>Intersection over union.</returns>
        public readonly float IoU(RectangleInt rectangle) => IoU(this, rectangle);

        /// <summary>
        /// Returns intersection area divided by union area.
        /// Degenerate or disjoint rectangles return zero.
        /// </summary>
        /// <param name="a">First rectangle.</param>
        /// <param name="b">Second rectangle.</param>
        /// <returns>Intersection over union.</returns>
        public static float IoU(RectangleInt a, RectangleInt b)
        {
            if (a.width <= 0 || a.height <= 0 || b.width <= 0 || b.height <= 0) return 0;
            double left = Math.Max((double)a.x, b.x);
            double top = Math.Max((double)a.y, b.y);
            double right = Math.Min((double)a.x + a.width, (double)b.x + b.width);
            double bottom = Math.Min((double)a.y + a.height, (double)b.y + b.height);
            double intersection = Math.Max(0, right - left) * Math.Max(0, bottom - top);
            if (intersection == 0) return 0;
            double union = (double)a.width * a.height + (double)b.width * b.height - intersection;
            return (float)(intersection / union);
        }

        /// <summary>
        /// Grows each dimension by its fractional increase around the center.
        /// Dimension increases and integer half-displacements are truncated toward zero.
        /// </summary>
        /// <param name="kx">Fractional increase in width.</param>
        /// <param name="ky">Fractional increase in height.</param>
        /// <returns>Resized rectangle.</returns>
        public readonly RectangleInt Scale(float kx = 0.0f, float ky = 0.0f)
        {
            int dx = (int)(width * kx);
            int dy = (int)(height * ky);
            return new RectangleInt(x - dx / 2, y - dy / 2, width + dx, height + dy);
        }

        /// <summary>
        /// Returns a centered square whose side is the original diagonal length.
        /// The diagonal and integer half-displacements are truncated toward zero.
        /// </summary>
        /// <returns>Square rectangle.</returns>
        public readonly RectangleInt Scale()
        {
            int side = (int)Math.Sqrt((double)width * width + (double)height * height);
            return new RectangleInt(x - (side - width) / 2, y - (side - height) / 2, side, side);
        }

        /// <summary>
        /// Normalizes this rectangle's negative dimensions and clips it to the bounds.
        /// Disjoint results retain the clipped origin and have nonnegative dimensions.
        /// </summary>
        /// <param name="bounds">Clipping rectangle, used without normalization.</param>
        /// <returns>Clipped rectangle.</returns>
        public readonly RectangleInt Clamp(RectangleInt bounds) => Clamp(this, bounds);

        /// <summary>
        /// Normalizes the first rectangle's negative dimensions and clips it to the second.
        /// </summary>
        /// <param name="first">Rectangle to normalize and clip.</param>
        /// <param name="second">Clipping rectangle, used without normalization.</param>
        /// <returns>Clipped rectangle.</returns>
        public static RectangleInt Clamp(RectangleInt first, RectangleInt second)
        {
            // Widen before adding or normalizing, including int.MinValue dimensions.
            long x2 = (long)first.x + first.width;
            long y2 = (long)first.y + first.height;
            long left = Math.Max(Math.Min(first.x, x2), second.x);
            long top = Math.Max(Math.Min(first.y, y2), second.y);
            long right = Math.Min(Math.Max(first.x, x2), (long)second.x + second.width);
            long bottom = Math.Min(Math.Max(first.y, y2), (long)second.y + second.height);
            return new RectangleInt((int)left, (int)top,
                (int)Math.Max(0, right - left), (int)Math.Max(0, bottom - top));
        }
        #endregion

        #region Conversions
        /// <summary>
        /// Rounds each floating-point component toward positive infinity.
        /// </summary>
        /// <param name="rectangle">Floating-point rectangle.</param>
        /// <returns>Integer rectangle.</returns>
        public static RectangleInt Ceiling(RectangleFloat rectangle)
        {
            return new RectangleInt(unchecked((int)Math.Ceiling(rectangle.X)),
                unchecked((int)Math.Ceiling(rectangle.Y)),
                unchecked((int)Math.Ceiling(rectangle.Width)),
                unchecked((int)Math.Ceiling(rectangle.Height)));
        }

        /// <summary>
        /// Rounds each floating-point component to the nearest integer, with midpoint ties to even.
        /// </summary>
        /// <param name="rectangle">Floating-point rectangle.</param>
        /// <returns>Integer rectangle.</returns>
        public static RectangleInt Round(RectangleFloat rectangle)
        {
            return new RectangleInt(unchecked((int)Math.Round(rectangle.X)),
                unchecked((int)Math.Round(rectangle.Y)),
                unchecked((int)Math.Round(rectangle.Width)),
                unchecked((int)Math.Round(rectangle.Height)));
        }

        /// <summary>
        /// Truncates each floating-point component toward zero.
        /// </summary>
        /// <param name="rectangle">Floating-point rectangle.</param>
        /// <returns>Integer rectangle.</returns>
        public static RectangleInt Truncate(RectangleFloat rectangle)
        {
            return new RectangleInt(unchecked((int)rectangle.X), unchecked((int)rectangle.Y),
                unchecked((int)rectangle.Width), unchecked((int)rectangle.Height));
        }

        /// <summary>
        /// Converts a System.Drawing.Rectangle without changing its components.
        /// </summary>
        /// <param name="rectangle">System.Drawing rectangle.</param>
        public static implicit operator RectangleInt(System.Drawing.Rectangle rectangle)
        {
            return new RectangleInt(rectangle.X, rectangle.Y, rectangle.Width, rectangle.Height);
        }

        /// <summary>
        /// Converts to System.Drawing.Rectangle without changing the components.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        public static implicit operator System.Drawing.Rectangle(RectangleInt rectangle)
        {
            return new System.Drawing.Rectangle(rectangle.x, rectangle.y, rectangle.width, rectangle.height);
        }
        #endregion

        #region Equality and overrides
        /// <summary>
        /// Compares all four components.
        /// </summary>
        /// <param name="a">First rectangle.</param>
        /// <param name="b">Second rectangle.</param>
        /// <returns>Boolean.</returns>
        public static bool operator ==(RectangleInt a, RectangleInt b)
        {
            return a.x == b.x && a.y == b.y && a.width == b.width && a.height == b.height;
        }

        /// <summary>
        /// Tests whether any component differs.
        /// </summary>
        /// <param name="a">First rectangle.</param>
        /// <param name="b">Second rectangle.</param>
        /// <returns>Boolean.</returns>
        public static bool operator !=(RectangleInt a, RectangleInt b) => !(a == b);

        /// <summary>
        /// Tests whether another rectangle has the same components.
        /// </summary>
        /// <param name="other">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public readonly bool Equals(RectangleInt other) => this == other;

        /// <summary>
        /// Tests whether an object is an equal RectangleInt.
        /// </summary>
        /// <param name="obj">Object.</param>
        /// <returns>Boolean.</returns>
        public override readonly bool Equals(object obj) => obj is RectangleInt other && Equals(other);

        /// <summary>
        /// Returns a hash code based on all four components.
        /// </summary>
        /// <returns>Hash code.</returns>
        public override readonly int GetHashCode()
        {
            unchecked
            {
                int hash = x.GetHashCode();
                hash = hash * 397 ^ y.GetHashCode();
                hash = hash * 397 ^ width.GetHashCode();
                return hash * 397 ^ height.GetHashCode();
            }
        }

        /// <summary>
        /// Returns the location and dimensions in System.Drawing format.
        /// </summary>
        /// <returns>Text representation.</returns>
        public override readonly string ToString()
        {
            return string.Format("{{X={0},Y={1},Width={2},Height={3}}}", x, y, width, height);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Returns a copy of this rectangle.
        /// </summary>
        /// <returns>Rectangle copy.</returns>
        public readonly RectangleInt Clone() => this;

        /// <summary>
        /// Returns a boxed copy of this rectangle.
        /// </summary>
        /// <returns>Rectangle copy.</returns>
        readonly object ICloneable.Clone() => Clone();
        #endregion
    }
}
