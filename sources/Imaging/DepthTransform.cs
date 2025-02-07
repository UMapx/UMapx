using System;
using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Uses for editing and transforming depth maps.
    /// </summary>
    public static class DepthTransform
    {
        #region Rotate
        /// <summary>
        /// Rotates depth by rotation value.
        /// </summary>
        /// <param name="depth">Matrix</param>
        /// <param name="rotation">Rotation</param>
        /// <returns>Matrix</returns>
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
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        private static ushort[,] Rotate90(ushort[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            ushort[,] H = new ushort[w, h];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[j, i] = input[h - i - 1, w - j - 1];
                }
            }

            return H;
        }
        /// <summary>
        /// Rotates the depth by 180 degrees.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        private static ushort[,] Rotate270(ushort[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            ushort[,] H = new ushort[w, h];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[j, i] = input[i, w - j - 1];
                }
            }

            return H;
        }
        #endregion
        
        /// <summary>
        /// Rotates depth by angle.
        /// </summary>
        /// <param name="depth">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <returns>Matrix</returns>
        public static ushort[,] Rotate(this ushort[,] depth, float angle)
        {
            return Rotate(depth, angle, 0);
        }
        /// <summary>
        /// Rotates depth by angle.
        /// </summary>
        /// <param name="depth">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="color">Background color</param>
        /// <returns>Matrix</returns>
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
            // coordinates of source points and cooefficiens
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
                            // get Y cooefficient
                            k1 = Kernel.Bicubic((float)(dy - n));

                            oy2 = oy1 + n;
                            if (oy2 < 0)
                                oy2 = 0;
                            if (oy2 > ymax)
                                oy2 = ymax;

                            for (int m = -1; m < 3; m++)
                            {
                                // get X cooefficient
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
        /// <param name="depth">Matrix</param>
        /// <param name="direction">Direction</param>
        /// <returns>Matrix</returns>
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
        /// <param name="depth">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// <param name="depth">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// <param name="depth">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// <param name="depth">Depth</param>
        /// <param name="rectangle">Rectangle</param>
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
        /// <param name="input">Matrix</param>
        /// <param name="size">Size</param>
        /// <returns>Matrix</returns>
        public static ushort[,] Resize(this ushort[,] input, Size size)
        {
            // get source size
            int width = input.GetLength(1);
            int height = input.GetLength(0);

            int w = size.Width;
            int h = size.Height;

            float xFactor = (float)width / w;
            float yFactor = (float)height / h;

            // coordinates of source points and cooefficiens
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
                        // get Y cooefficient
                        k1 = Kernel.Bicubic((float)(dy - n));

                        oy2 = oy1 + n;
                        if (oy2 < 0)
                            oy2 = 0;
                        if (oy2 > ymax)
                            oy2 = ymax;

                        for (int m = -1; m < 3; m++)
                        {
                            // get X cooefficient
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
        /// <param name="a">Matrix</param>
        /// <param name="h">The number of positions to which a shift in height occurs</param>
        /// <param name="w">The number of positions by which the shift occurs in width</param>
        /// <returns>Matrix</returns>
        public static ushort[,] Shift(this ushort[,] a, int w, int h)
        {
            int l0 = a.GetLength(0), l1 = a.GetLength(1);
            ushort[,] temp = new ushort[l0, l1];
            int i, j;

            for (i = 0; i < l0; i++)
            {
                for (j = 0; j < l1; j++)
                {
                    temp[i, j] = a[Maths.Mod(i - h, l1), Maths.Mod(j - w, l0)];
                }
            }
            return temp;
        }
        #endregion

        #region Merge
        /// <summary>
        /// Merges two depths.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="b">Matrix</param>
        public static void Merge(this ushort[,] a, ushort[,] b)
        {
            var rectangle = new Rectangle(0, 0, b.GetLength(1), b.GetLength(0));
            Merge(a, b, rectangle);
        }
        /// <summary>
        /// Merges two depths.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="b">Matrix</param>
        /// <param name="rectangle">Rectangle</param>
        public static void Merge(this ushort[,] a, ushort[,] b, Rectangle rectangle)
        {
            ushort[,] c = Resize(b, 
                new Size(rectangle.Width, rectangle.Height));

            int h = Math.Min(rectangle.Height, a.GetLength(0) - rectangle.Y);
            int w = Math.Min(rectangle.Width,  a.GetLength(1) - rectangle.X);

            for (int i = rectangle.Y; i < h; i++)
            {
                for (int j = rectangle.X; j < w; j++)
                {
                    a[i, j] = c[i - rectangle.Y, j - rectangle.X];
                }
            }
        }
        #endregion

        #region Equalize
        /// <summary>
        /// Equalizes histogram of the depth.
        /// </summary>
        /// <param name="depth">Depth</param>
        /// <returns>Matrix</returns>
        public static ushort[,] Equalize(this ushort[,] depth)
        {
            var width = depth.GetLength(1);
            var height = depth.GetLength(0);
            var hist = ushort.MaxValue + 1;

            // histogram
            var H = new ushort[hist];

            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    H[depth[y, x]]++;
                }
            }

            // cdf
            var factor = ushort.MaxValue / (float)(height * width);
            var cdf = new float[hist];

            // recursion
            cdf[0] = H[0];

            for (int i = 1; i < hist; i++)
            {
                cdf[i] = H[i] + cdf[i - 1];
            }

            // equalization
            var output = new ushort[height, width];

            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    output[y, x] = (ushort)(cdf[depth[y, x]] * factor);
                }
            }

            return output;
        }
        #endregion
    }
}
