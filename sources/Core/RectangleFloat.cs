using System;
using System.Numerics;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a rectangle with single-precision floating-point coordinates and dimensions.
    /// </summary>
    /// <remarks>
    /// Standard geometry members follow System.Drawing.RectangleF semantics.
    /// Additional arithmetic members return new values without changing this rectangle.
    /// </remarks>
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
        public float X { readonly get => x; set => x = value; }

        /// <summary>
        /// Gets or sets the top coordinate.
        /// </summary>
        public float Y { readonly get => y; set => y = value; }

        /// <summary>
        /// Gets or sets the width.
        /// </summary>
        public float Width { readonly get => width; set => width = value; }

        /// <summary>
        /// Gets or sets the height.
        /// </summary>
        public float Height { readonly get => height; set => height = value; }

        /// <summary>
        /// Gets or sets the upper-left corner.
        /// </summary>
        public PointFloat Location
        {
            readonly get => new PointFloat(x, y);
            set { x = value.X; y = value.Y; }
        }

        /// <summary>
        /// Gets or sets the dimensions.
        /// </summary>
        public SizeFloat Size
        {
            readonly get => new SizeFloat(width, height);
            set { width = value.Width; height = value.Height; }
        }

        /// <summary>
        /// Gets the left coordinate.
        /// </summary>
        public readonly float Left => x;

        /// <summary>
        /// Gets the top coordinate.
        /// </summary>
        public readonly float Top => y;

        /// <summary>
        /// Gets the right coordinate (X + Width).
        /// </summary>
        public readonly float Right => x + width;

        /// <summary>
        /// Gets the bottom coordinate (Y + Height).
        /// </summary>
        public readonly float Bottom => y + height;

        /// <summary>
        /// Tests whether either dimension is zero or negative, regardless of location.
        /// </summary>
        /// <remarks>
        /// Uses the same rule as RectangleInt. Compare with Empty to test whether all components are zero.
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
        public readonly bool Contains(float x, float y)
        {
            return this.x <= x && x < Right && this.y <= y && y < Bottom;
        }

        /// <summary>
        /// Tests whether a point lies inside, excluding the right and bottom edges.
        /// </summary>
        /// <param name="point">Point.</param>
        /// <returns>Boolean.</returns>
        public readonly bool Contains(PointFloat point) => Contains(point.X, point.Y);

        /// <summary>
        /// Tests whether the specified rectangle is entirely contained within this rectangle.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public readonly bool Contains(RectangleFloat rectangle)
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
        public readonly bool IntersectsWith(RectangleFloat rectangle)
        {
            return rectangle.x < Right && x < rectangle.Right
                && rectangle.y < Bottom && y < rectangle.Bottom;
        }

        /// <summary>
        /// Replaces this rectangle with its intersection with the specified rectangle.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        public void Intersect(RectangleFloat rectangle)
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
        /// <param name="rectangle">Rectangle.</param>
        /// <param name="x">Horizontal amount per side.</param>
        /// <param name="y">Vertical amount per side.</param>
        /// <returns>Expanded rectangle.</returns>
        public static RectangleFloat Inflate(RectangleFloat rectangle, float x, float y)
        {
            rectangle.Inflate(x, y);
            return rectangle;
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
        /// <param name="point">Displacement.</param>
        public void Offset(PointFloat point) => Offset(point.X, point.Y);
        #endregion

        #region Arithmetic
        /// <summary>
        /// Returns a translated copy of this rectangle.
        /// </summary>
        /// <param name="point">Displacement to add.</param>
        /// <returns>Translated rectangle.</returns>
        public readonly RectangleFloat Add(PointFloat point)
        {
            return new RectangleFloat(x + point.X, y + point.Y, width, height);
        }

        /// <summary>
        /// Returns a copy translated by the negative of the displacement.
        /// </summary>
        /// <param name="point">Displacement to subtract.</param>
        /// <returns>Translated rectangle.</returns>
        public readonly RectangleFloat Sub(PointFloat point)
        {
            return new RectangleFloat(x - point.X, y - point.Y, width, height);
        }

        /// <summary>
        /// Adds a displacement without changing the original rectangle.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        /// <param name="point">Displacement.</param>
        /// <returns>Translated rectangle.</returns>
        public static RectangleFloat operator +(RectangleFloat rectangle, PointFloat point) => rectangle.Add(point);

        /// <summary>
        /// Subtracts a displacement without changing the original rectangle.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        /// <param name="point">Displacement.</param>
        /// <returns>Translated rectangle.</returns>
        public static RectangleFloat operator -(RectangleFloat rectangle, PointFloat point) => rectangle.Sub(point);

        /// <summary>
        /// Returns translated copies of an array of rectangles.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <param name="point">Displacement to add.</param>
        /// <returns>New array of translated rectangles.</returns>
        public static RectangleFloat[] Add(RectangleFloat[] rectangles, PointFloat point)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = new RectangleFloat[rectangles.Length];
            for (int i = 0; i < result.Length; i++) result[i] = rectangles[i].Add(point);
            return result;
        }

        /// <summary>
        /// Returns copies of an array translated by the negative displacement.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <param name="point">Displacement to subtract.</param>
        /// <returns>New array of translated rectangles.</returns>
        public static RectangleFloat[] Sub(RectangleFloat[] rectangles, PointFloat point)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = new RectangleFloat[rectangles.Length];
            for (int i = 0; i < result.Length; i++) result[i] = rectangles[i].Sub(point);
            return result;
        }

        /// <summary>
        /// Returns corners clockwise: top-left, top-right, bottom-right, bottom-left.
        /// </summary>
        /// <returns>Four corner points.</returns>
        public readonly PointFloat[] ToPoints()
        {
            return new[] { new PointFloat(Left, Top), new PointFloat(Right, Top),
                new PointFloat(Right, Bottom), new PointFloat(Left, Bottom) };
        }

        /// <summary>
        /// Creates a rectangle from four ordered corners, using corners zero and two.
        /// </summary>
        /// <param name="points">Four corners in the order returned by ToPoints.</param>
        /// <returns>Rectangle.</returns>
        /// <exception cref="ArgumentNullException">The points array is null.</exception>
        /// <exception cref="ArgumentException">The array does not contain four points.</exception>
        public static RectangleFloat FromPoints(PointFloat[] points)
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
        public readonly PointFloat GetPoint() => Location;

        /// <summary>
        /// Returns Width multiplied by Height.
        /// </summary>
        /// <returns>Signed area.</returns>
        public readonly float Area() => width * height;

        /// <summary>
        /// Returns the product of the specified dimensions.
        /// </summary>
        /// <param name="size">Dimensions.</param>
        /// <returns>Signed area.</returns>
        public static float Area(SizeFloat size) => size.Width * size.Height;

        /// <summary>
        /// Returns the first non-empty rectangle with the largest signed area.
        /// Products are compared in double precision to avoid area overflow.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <returns>Selected rectangle, or Empty if no non-empty rectangle exists.</returns>
        public static RectangleFloat Max(params RectangleFloat[] rectangles)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = Empty;
            double bestArea = 0;
            bool found = false;
            foreach (var rectangle in rectangles)
            {
                if (rectangle.IsEmpty) continue;
                double area = (double)rectangle.width * rectangle.height;
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
        /// Products are compared in double precision to avoid area overflow.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <returns>Selected rectangle, or Empty if no non-empty rectangle exists.</returns>
        public static RectangleFloat Min(params RectangleFloat[] rectangles)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = Empty;
            double bestArea = 0;
            bool found = false;
            foreach (var rectangle in rectangles)
            {
                if (rectangle.IsEmpty) continue;
                double area = (double)rectangle.width * rectangle.height;
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
        /// The resulting side is the larger original dimension.
        /// </summary>
        /// <returns>Square rectangle.</returns>
        public readonly RectangleFloat ToBox()
        {
            float side = Math.Max(width, height);
            return new RectangleFloat(x - (side - width) / 2, y - (side - height) / 2, side, side);
        }

        /// <summary>
        /// Grows each dimension by the specified fraction around the center.
        /// A scale of zero preserves the rectangle; one doubles its dimensions.
        /// This overload preserves the aspect ratio instead of forming a square.
        /// </summary>
        /// <param name="scale">Fractional increase in both dimensions.</param>
        /// <returns>Resized rectangle.</returns>
        public readonly RectangleFloat ToBox(float scale)
        {
            float dx = width * scale;
            float dy = height * scale;
            return new RectangleFloat(x - dx / 2, y - dy / 2, width + dx, height + dy);
        }

        /// <summary>
        /// Returns square copies of the specified rectangles.
        /// </summary>
        /// <param name="rectangles">Rectangles.</param>
        /// <returns>New array of square rectangles.</returns>
        public static RectangleFloat[] ToBox(params RectangleFloat[] rectangles)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = new RectangleFloat[rectangles.Length];
            for (int i = 0; i < result.Length; i++) result[i] = rectangles[i].ToBox();
            return result;
        }

        /// <summary>
        /// Returns copies grown by the specified fraction around their centers.
        /// </summary>
        /// <param name="factor">Fractional increase in both dimensions.</param>
        /// <param name="rectangles">Rectangles.</param>
        /// <returns>New array of resized rectangles.</returns>
        public static RectangleFloat[] ToBox(float factor, params RectangleFloat[] rectangles)
        {
            if (rectangles == null) throw new ArgumentNullException(nameof(rectangles));
            var result = new RectangleFloat[rectangles.Length];
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
        public readonly float IoU(RectangleFloat rectangle) => IoU(this, rectangle);

        /// <summary>
        /// Returns intersection area divided by union area.
        /// Degenerate or disjoint rectangles return zero.
        /// </summary>
        /// <param name="a">First rectangle.</param>
        /// <param name="b">Second rectangle.</param>
        /// <returns>Intersection over union.</returns>
        public static float IoU(RectangleFloat a, RectangleFloat b)
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
        /// Zero preserves that dimension; one doubles it.
        /// </summary>
        /// <param name="kx">Fractional increase in width.</param>
        /// <param name="ky">Fractional increase in height.</param>
        /// <returns>Resized rectangle.</returns>
        public readonly RectangleFloat Scale(float kx = 0.0f, float ky = 0.0f)
        {
            float dx = width * kx;
            float dy = height * ky;
            return new RectangleFloat(x - dx / 2, y - dy / 2, width + dx, height + dy);
        }

        /// <summary>
        /// Returns a centered square whose side is the original diagonal length.
        /// Intermediate products use double precision.
        /// </summary>
        /// <returns>Square rectangle.</returns>
        public readonly RectangleFloat Scale()
        {
            float side = (float)Math.Sqrt((double)width * width + (double)height * height);
            return new RectangleFloat(x - (side - width) / 2, y - (side - height) / 2, side, side);
        }

        /// <summary>
        /// Normalizes this rectangle's negative dimensions and clips it to the bounds.
        /// Disjoint results retain the clipped origin and have nonnegative dimensions.
        /// </summary>
        /// <param name="bounds">Clipping rectangle, used without normalization.</param>
        /// <returns>Clipped rectangle.</returns>
        public readonly RectangleFloat Clamp(RectangleFloat bounds) => Clamp(this, bounds);

        /// <summary>
        /// Normalizes the first rectangle's negative dimensions and clips it to the second.
        /// </summary>
        /// <param name="first">Rectangle to normalize and clip.</param>
        /// <param name="second">Clipping rectangle, used without normalization.</param>
        /// <returns>Clipped rectangle.</returns>
        public static RectangleFloat Clamp(RectangleFloat first, RectangleFloat second)
        {
            if (first.width < 0) { first.x += first.width; first.width = -first.width; }
            if (first.height < 0) { first.y += first.height; first.height = -first.height; }
            float left = Math.Max(first.Left, second.Left);
            float top = Math.Max(first.Top, second.Top);
            float right = Math.Min(first.Right, second.Right);
            float bottom = Math.Min(first.Bottom, second.Bottom);
            return new RectangleFloat(left, top, Math.Max(0, right - left), Math.Max(0, bottom - top));
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
        public readonly Vector4 ToVector4() => new Vector4(x, y, width, height);

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
        /// <param name="rectangle">Integer rectangle.</param>
        public static implicit operator RectangleFloat(RectangleInt rectangle)
        {
            return new RectangleFloat(rectangle.X, rectangle.Y, rectangle.Width, rectangle.Height);
        }

        /// <summary>
        /// Converts a System.Drawing.RectangleF without changing its components.
        /// </summary>
        /// <param name="rectangle">System.Drawing rectangle.</param>
        public static implicit operator RectangleFloat(System.Drawing.RectangleF rectangle)
        {
            return new RectangleFloat(rectangle.X, rectangle.Y, rectangle.Width, rectangle.Height);
        }

        /// <summary>
        /// Converts to System.Drawing.RectangleF without changing the components.
        /// </summary>
        /// <param name="rectangle">Rectangle.</param>
        public static implicit operator System.Drawing.RectangleF(RectangleFloat rectangle)
        {
            return new System.Drawing.RectangleF(rectangle.x, rectangle.y, rectangle.width, rectangle.height);
        }
        #endregion

        #region Equality and overrides
        /// <summary>
        /// Compares all four components.
        /// </summary>
        /// <param name="a">First rectangle.</param>
        /// <param name="b">Second rectangle.</param>
        /// <returns>Boolean.</returns>
        public static bool operator ==(RectangleFloat a, RectangleFloat b)
        {
            return a.x == b.x && a.y == b.y && a.width == b.width && a.height == b.height;
        }

        /// <summary>
        /// Tests whether any component differs.
        /// </summary>
        /// <param name="a">First rectangle.</param>
        /// <param name="b">Second rectangle.</param>
        /// <returns>Boolean.</returns>
        public static bool operator !=(RectangleFloat a, RectangleFloat b) => !(a == b);

        /// <summary>
        /// Tests whether another rectangle has the same components.
        /// </summary>
        /// <param name="other">Rectangle.</param>
        /// <returns>Boolean.</returns>
        public readonly bool Equals(RectangleFloat other) => this == other;

        /// <summary>
        /// Tests whether an object is an equal RectangleFloat.
        /// </summary>
        /// <param name="obj">Object.</param>
        /// <returns>Boolean.</returns>
        public override readonly bool Equals(object obj) => obj is RectangleFloat other && Equals(other);

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
        public readonly RectangleFloat Clone() => this;

        /// <summary>
        /// Returns a boxed copy of this rectangle.
        /// </summary>
        /// <returns>Rectangle copy.</returns>
        readonly object ICloneable.Clone() => Clone();
        #endregion
    }
}
