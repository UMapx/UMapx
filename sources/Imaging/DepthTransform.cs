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
        #region Public methods

        /// <summary>
        /// Converts Bitmap into unshort matrix.
        /// </summary>
        /// <param name="depth">Image</param>
        /// <returns>Bitmap</returns>
        public unsafe static ushort[,] FromBitmap(this Bitmap depth)
        {
            var bmData = depth.Lock32bpp();
            var width = bmData.Width;
            var height = bmData.Height;
            var stride = bmData.Stride;
            var src = (byte*)bmData.Scan0.ToPointer();
            var output = new ushort[height, width];
            
            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    var k = x * 3 + y * stride;
                    output[y, x] = (ushort)RGB.Average(src[k + 0], src[k + 1], src[k + 2]);
                }
            }

            depth.Unlock(bmData);
            return output;
        }
        
        /// <summary>
        /// Converts ushort matrix into Bitmap.
        /// </summary>
        /// <param name="depth">Matrix</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap ToBitmap(this ushort[,] depth)
        {
            var width = depth.GetLength(1);
            var height = depth.GetLength(0);
            var rectangle = new Rectangle(0, 0, width, height);
            var bitmap = new Bitmap(width, height);
            var bmData = bitmap.LockBits(rectangle, ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            var dst = (byte*)bmData.Scan0.ToPointer();
            var stride = bmData.Stride;

            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    var k = x * 3 + y * stride;
                    dst[k + 0] = dst[k + 1] = dst[k + 2] = Maths.Byte((float)depth[y, x] / byte.MaxValue);
                }
            }

            bitmap.Unlock(bmData);
            return bitmap;
        }

        /// <summary>
        /// Converts the depth to the matrix.
        /// </summary>
        /// <param name="depth">Depth</param>
        /// <returns>Matrix</returns>
        public static float[,] ToFloat(this ushort[,] depth)
        {
            int h = depth.GetLength(0);
            int w = depth.GetLength(1);
            float[,] output = new float[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    output[i, j] = depth[i, j] / (float)ushort.MaxValue;
                }
            }

            return output;
        }

        /// <summary>
        /// Converts the ushort matrix to the float matrix.
        /// </summary>
        /// <param name="depth">Matrix</param>
        /// <returns>Matrix</returns>
        public static ushort[,] FromFloat(this float[,] depth)
        {
            int h = depth.GetLength(0);
            int w = depth.GetLength(1);
            ushort[,] output = new ushort[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    output[i, j] = (ushort)(depth[i, j] * ushort.MaxValue);
                }
            }

            return output;
        }

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

        /// <summary>
        /// Returns resized matrix.
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
    }
}
