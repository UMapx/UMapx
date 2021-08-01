using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Uses for editing and transforming depth maps.
    /// </summary>
    public static class DepthTransform
    {
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
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <returns>Matrix</returns>
        public static ushort[,] Resize(this ushort[,] input, int w, int h)
        {
            // get source size
            int width = input.GetLength(1);
            int height = input.GetLength(0);

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

                    H[y, x] = (ushort)g;
                }
            }

            return H;
        }
        #endregion

        #region Canvas
        /// <summary>
        /// Rotates the depth by 90 degrees.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        public static ushort[,] Rotate90(this ushort[,] input)
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
        public static ushort[,] Rotate180(this ushort[,] input)
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
        public static ushort[,] Rotate270(this ushort[,] input)
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
    }
}
