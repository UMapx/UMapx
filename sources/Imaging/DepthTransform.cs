using System;
using System.Drawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Used to edit and transform depth maps.
    /// </summary>
    public static partial class DepthTransform
    {
        #region Rotate
        /// <summary>
        /// Rotates depth by rotation value.
        /// </summary>
        /// <param name="depth">Matrix.</param>
        /// <param name="rotation">Rotation.</param>
        /// <returns>Matrix.</returns>
        public static ushort[,] Rotate(this ushort[,] depth, RotationMode rotation)
        {
            switch (rotation)
            {
                case RotationMode.R0:
                    return depth;
                case RotationMode.R90:
                    return Rotate90(depth);
                case RotationMode.R180:
                    return Rotate180(depth);
                case RotationMode.R270:
                    return Rotate270(depth);
                default:
                    return depth;
            }
        }

        #region Private
        /// <summary>
        /// Rotates the depth by 90 degrees.
        /// </summary>
        /// <param name="input">Matrix.</param>
        /// <returns>Matrix.</returns>
        private static ushort[,] Rotate90(ushort[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            ushort[,] H = new ushort[w, h];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[j, h - 1 - i] = input[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Rotates the depth by 180 degrees.
        /// </summary>
        /// <param name="input">Matrix.</param>
        /// <returns>Matrix.</returns>
        private static ushort[,] Rotate180(ushort[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            ushort[,] H = new ushort[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[i, j] = input[h - i - 1, w - j - 1];
                }
            }

            return H;
        }
        /// <summary>
        /// Rotates the depth by 270 degrees.
        /// </summary>
        /// <param name="input">Matrix.</param>
        /// <returns>Matrix.</returns>
        private static ushort[,] Rotate270(ushort[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            ushort[,] H = new ushort[w, h];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[w - 1 - j, i] = input[i, j];
                }
            }

            return H;
        }
        #endregion
        
        /// <summary>
        /// Rotates depth by angle.
        /// </summary>
        /// <param name="depth">Matrix.</param>
        /// <param name="angle">Angle.</param>
        /// <returns>Matrix.</returns>
        public static ushort[,] Rotate(this ushort[,] depth, float angle)
        {
            return Rotate(depth, angle, 0);
        }
        /// <summary>
        /// Rotates depth by angle.
        /// </summary>
        /// <param name="depth">Matrix.</param>
        /// <param name="angle">Angle.</param>
        /// <param name="color">Background color.</param>
        /// <returns>Matrix.</returns>
        public static ushort[,] Rotate(this ushort[,] depth, float angle, ushort color)
        {
            // get source image size
            int width = depth.GetLength(1);
            int height = depth.GetLength(0);
            float oldXradius = (float)(width - 1) / 2;
            float oldYradius = (float)(height - 1) / 2;

            // get destination image size
            int newWidth = width;
            int newHeight = height;
            float newXradius = (float)(newWidth - 1) / 2;
            float newYradius = (float)(newHeight - 1) / 2;

            // angle's sine and cosine
            float angleRad = -angle * Maths.Pi / 180.0f;
            float angleCos = Maths.Cos(angleRad);
            float angleSin = Maths.Sin(angleRad);

            // destination pixel's coordinate relative to image center
            float cx, cy;
            // coordinates of source points and coefficients
            float ox, oy, dx, dy, k1, k2;
            int ox1, oy1, ox2, oy2;
            // destination pixel values
            float g;
            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;
            // output
            ushort[,] H = new ushort[newHeight, newWidth];

            // grayscale
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                cx = -newXradius;

                for (int x = 0; x < newWidth; x++)
                {
                    // coordinates of source point
                    ox = angleCos * cx + angleSin * cy + oldXradius;
                    oy = -angleSin * cx + angleCos * cy + oldYradius;

                    ox1 = (int)ox;
                    oy1 = (int)oy;

                    // validate source pixel's coordinates
                    if ((ox1 < 0) || (oy1 < 0) || (ox1 >= width) || (oy1 >= height))
                    {
                        // fill destination image with filler
                        H[y, x] = color;
                    }
                    else
                    {
                        dx = ox - ox1;
                        dy = oy - oy1;

                        // initial pixel value
                        g = 0;

                        for (int n = -1; n < 3; n++)
                        {
                            // get Y coefficient
                            k1 = Kernel.Bicubic((float)(dy - n));

                            oy2 = oy1 + n;
                            if (oy2 < 0)
                                oy2 = 0;
                            if (oy2 > ymax)
                                oy2 = ymax;

                            for (int m = -1; m < 3; m++)
                            {
                                // get X coefficient
                                k2 = k1 * Kernel.Bicubic((float)(m - dx));

                                ox2 = ox1 + m;
                                if (ox2 < 0)
                                    ox2 = 0;
                                if (ox2 > xmax)
                                    ox2 = xmax;

                                g += k2 * depth[oy2, ox2];
                            }
                        }
                        H[y, x] = (ushort)Maths.Range(g, ushort.MinValue, ushort.MaxValue);
                    }
                    cx++;
                }
                cy++;
            }

            return H;
        }
        #endregion

        #region Flip
        /// <summary>
        /// Flips depth by direction.
        /// </summary>
        /// <param name="depth">Matrix.</param>
        /// <param name="direction">Direction.</param>
        /// <returns>Matrix.</returns>
        public static ushort[,] Flip(this ushort[,] depth, Direction direction)
        {
            switch (direction)
            {
                case Direction.Horizontal:
                    return FlipX(depth);
                case Direction.Vertical:
                    return FlipY(depth);
                case Direction.Both:
                    return FlipXY(depth);
                default:
                    return depth;
            }
        }

        #region Private
        /// <summary>
        /// Flips depth by X axis.
        /// </summary>
        /// <param name="depth">Matrix.</param>
        /// <returns>Matrix.</returns>
        private static ushort[,] FlipX(ushort[,] depth)
        {
            int h = depth.GetLength(0);
            int w = depth.GetLength(1);

            ushort[,] H = new ushort[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[i, j] = depth[i, w - j - 1];
                }
            }

            return H;
        }
        /// <summary>
        /// Flips depth by Y axis.
        /// </summary>
        /// <param name="depth">Matrix.</param>
        /// <returns>Matrix.</returns>
        private static ushort[,] FlipY(ushort[,] depth)
        {
            int h = depth.GetLength(0);
            int w = depth.GetLength(1);

            ushort[,] H = new ushort[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[i, j] = depth[h - i - 1, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Flips depth by XY axis.
        /// </summary>
        /// <param name="depth">Matrix.</param>
        /// <returns>Matrix.</returns>
        private static ushort[,] FlipXY(ushort[,] depth)
        {
            int h = depth.GetLength(0);
            int w = depth.GetLength(1);

            ushort[,] H = new ushort[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[i, j] = depth[h - i - 1, w - j - 1];
                }
            }

            return H;
        }
        #endregion
        #endregion

        #region Crop
        /// <summary>
        /// Crops the depth.
        /// </summary>
        /// <param name="depth">Depth.</param>
        /// <param name="rectangle">Rectangle.</param>
        /// <returns></returns>
        public static ushort[,] Crop(this ushort[,] depth, Rectangle rectangle)
        {
            // image params
            int width = depth.GetLength(1);
            int height = depth.GetLength(0);

            // check section params
            int x = Maths.Range(rectangle.X, 0, width);
            int y = Maths.Range(rectangle.Y, 0, height);
            int w = Maths.Range(rectangle.Width, 0, width - x);
            int h = Maths.Range(rectangle.Height, 0, height - y);

            // exception
            if (x == 0 &&
                y == 0 &&
                w == 0 &&
                h == 0) return depth;

            // output
            var output = new ushort[h, w];

            for (int i = 0; i < w; i++)
            {
                for (int j = 0; j < h; j++)
                {
                    output[j, i] = depth[y + j, x + i];
                }
            }

            return output;
        }
        #endregion

        #region Resize
        /// <summary>
        /// Resizes the depth.
        /// </summary>
        /// <param name="input">Matrix.</param>
        /// <param name="size">Size.</param>
        /// <returns>Matrix.</returns>
        public static ushort[,] Resize(this ushort[,] input, Size size)
        {
            // get source size
            int width = input.GetLength(1);
            int height = input.GetLength(0);

            int w = size.Width;
            int h = size.Height;

            float xFactor = (float)width / w;
            float yFactor = (float)height / h;

            // coordinates of source points and coefficients
            float ox, oy, dx, dy, k1, k2;
            int ox1, oy1, ox2, oy2;
            float g;

            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;

            // output
            ushort[,] H = new ushort[h, w];

            // grayscale
            for (int y = 0; y < h; y++)
            {
                // Y coordinates
                oy = y * yFactor - 0.5f;
                oy1 = (int)oy;
                dy = oy - oy1;

                for (int x = 0; x < w; x++)
                {
                    // X coordinates
                    ox = x * xFactor - 0.5f;
                    ox1 = (int)ox;
                    dx = ox - ox1;

                    // initial pixel value
                    g = 0;

                    for (int n = -1; n < 3; n++)
                    {
                        // get Y coefficient
                        k1 = Kernel.Bicubic((float)(dy - n));

                        oy2 = oy1 + n;
                        if (oy2 < 0)
                            oy2 = 0;
                        if (oy2 > ymax)
                            oy2 = ymax;

                        for (int m = -1; m < 3; m++)
                        {
                            // get X coefficient
                            k2 = k1 * Kernel.Bicubic((float)(m - dx));

                            ox2 = ox1 + m;
                            if (ox2 < 0)
                                ox2 = 0;
                            if (ox2 > xmax)
                                ox2 = xmax;

                            g += k2 * input[oy2, ox2];
                        }
                    }

                    H[y, x] = (ushort)Maths.Range(g, ushort.MinValue, ushort.MaxValue);
                }
            }

            return H;
        }
        #endregion

        #region Shift
        /// <summary>
        /// Shifts the depth.
        /// </summary>
        /// <param name="a">Matrix.</param>
        /// <param name="h">The number of positions to which a shift in height occurs.</param>
        /// <param name="w">The number of positions by which the shift occurs in width.</param>
        /// <returns>Matrix.</returns>
        public static ushort[,] Shift(this ushort[,] a, int w, int h)
        {
            int l0 = a.GetLength(0), l1 = a.GetLength(1);
            ushort[,] temp = new ushort[l0, l1];
            int i, j;

            for (i = 0; i < l0; i++)
            {
                for (j = 0; j < l1; j++)
                {
                    temp[i, j] = a[Maths.Mod(i - h, l0), Maths.Mod(j - w, l1)];
                }
            }
            return temp;
        }
        #endregion

        #region Merge
        /// <summary>
        /// Copies a depth map into the top-left corner of another depth map, clipping at its edges.
        /// </summary>
        /// <param name="a">Matrix.</param>
        /// <param name="b">Matrix.</param>
        public static void Merge(this ushort[,] a, ushort[,] b)
        {
            var rectangle = new Rectangle(0, 0, b.GetLength(1), b.GetLength(0));
            Merge(a, b, rectangle);
        }
        /// <summary>
        /// Resizes a depth map to a placement rectangle and copies its visible part into another map.
        /// </summary>
        /// <param name="a">Destination map, modified in place.</param>
        /// <param name="b">Source map. Equal source and placement sizes preserve samples exactly.</param>
        /// <param name="rectangle">Placement in column/row coordinates, with nonnegative dimensions.</param>
        /// <remarks>Clipping retains the source coordinates relative to the full placement.
        /// Empty or nonintersecting placements are no-ops. Self-merges read the original samples.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">A placement dimension is negative.</exception>
        /// <exception cref="ArgumentException">A visible, nonempty placement has an empty source.</exception>
        public static void Merge(this ushort[,] a, ushort[,] b, Rectangle rectangle)
        {
            if (rectangle.Width < 0 || rectangle.Height < 0)
                throw new ArgumentOutOfRangeException(nameof(rectangle), "Placement dimensions must be nonnegative.");

            int top = Math.Max(0, rectangle.Y), left = Math.Max(0, rectangle.X);
            // Widen before adding coordinates, even for placements entirely outside the map.
            long bottom = Math.Min(a.GetLength(0), (long)rectangle.Y + rectangle.Height);
            long right = Math.Min(a.GetLength(1), (long)rectangle.X + rectangle.Width);
            if (top >= bottom || left >= right) return;
            if (b.Length == 0)
                throw new ArgumentException("A nonempty placement requires source samples.", nameof(b));

            ushort[,] c = b.GetLength(1) == rectangle.Width && b.GetLength(0) == rectangle.Height
                ? (ReferenceEquals(a, b) ? (ushort[,])b.Clone() : b)
                : Resize(b, rectangle.Size);

            for (int i = top; i < bottom; i++)
            {
                for (int j = left; j < right; j++)
                {
                    a[i, j] = c[i - rectangle.Y, j - rectangle.X];
                }
            }
        }
        #endregion

        #region Equalize
        /// <summary>
        /// Equalizes a depth map using its inclusive cumulative histogram.
        /// </summary>
        /// <param name="depth">Unsigned 16-bit depth samples.</param>
        /// <returns>A map of the same shape, with each value mapped to
        /// floor(65535 * count(samples less than or equal to value) / population).
        /// Constant nonempty maps become 65535; empty maps retain their shape.</returns>
        public static ushort[,] Equalize(this ushort[,] depth)
        {
            var width = depth.GetLength(1);
            var height = depth.GetLength(0);
            var hist = ushort.MaxValue + 1;
            var output = new ushort[height, width];
            long population = depth.LongLength;
            if (population == 0) return output;

            // histogram
            var H = new long[hist];

            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    H[depth[y, x]]++;
                }
            }

            // Integer CDF ranks avoid both 16-bit count wraparound and float rounding
            // near output bin boundaries. The product fits Int64 for managed-array sizes.
            var lookup = new ushort[hist];
            long cumulative = 0;
            for (int i = 0; i < hist; i++)
            {
                cumulative += H[i];
                lookup[i] = (ushort)(cumulative * ushort.MaxValue / population);
            }

            // equalization
            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    output[y, x] = lookup[depth[y, x]];
                }
            }

            return output;
        }
        #endregion
    }
}
