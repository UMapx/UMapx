using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines a Canny edge detector.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Canny_edge_detector
    /// </remarks>
    /// </summary>
    [Serializable]
    public class CannyEdgeDetector : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private float lowThreshold;
        private float highThreshold;
        private int gaussianRadius;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the Canny edge detector.
        /// </summary>
        /// <param name="lowThreshold">Low threshold</param>
        /// <param name="highThreshold">High threshold</param>
        /// <param name="radius">Gaussian blur radius</param>
        public CannyEdgeDetector(float lowThreshold = 20f, float highThreshold = 60f, int radius = 2)
        {
            LowThreshold = lowThreshold;
            HighThreshold = highThreshold;
            Radius = radius;
        }
        /// <summary>
        /// Gets or sets the low threshold.
        /// </summary>
        public float LowThreshold
        {
            get
            {
                return this.lowThreshold;
            }
            set
            {
                this.lowThreshold = value;
            }
        }
        /// <summary>
        /// Gets or sets the high threshold.
        /// </summary>
        public float HighThreshold
        {
            get
            {
                return this.highThreshold;
            }
            set
            {
                this.highThreshold = value;
            }
        }
        /// <summary>
        /// Gets or sets the Gaussian blur radius.
        /// </summary>
        public int Radius
        {
            get
            {
                return this.gaussianRadius;
            }
            set
            {
                this.gaussianRadius = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            int width = bmSrc.Width;
            int height = bmSrc.Height;
            int stride = bmSrc.Stride;

            byte* src = (byte*)bmSrc.Scan0;
            byte* dst = (byte*)bmData.Scan0;

            float[,] gray = new float[height, width];
            float[,] blurred = new float[height, width];
            float[,] gradient = new float[height, width];
            float[,] direction = new float[height, width];

            // Convert to grayscale
            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    byte* p = &src[y * stride + x * 4];
                    float grayVal = 0.299f * p[2] + 0.587f * p[1] + 0.114f * p[0];
                    gray[y, x] = grayVal;
                }
            });

            // Apply Gaussian blur
            GaussianBlur(gray, blurred, gaussianRadius);

            // Sobel gradients
            float[,] gx = new float[height, width];
            float[,] gy = new float[height, width];

            // Sobel filter
            ApplySobel(blurred, gx, gy, gradient, direction);

            // Non-maximum suppression
            float[,] nms = new float[height, width];
            NonMaximumSuppression(gradient, direction, nms);

            // Hysteresis thresholding
            bool[,] edges = new bool[height, width];
            Hysteresis(nms, edges, lowThreshold, highThreshold);

            // Write output
            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    byte value = edges[y, x] ? (byte)255 : (byte)0;
                    byte* p = &dst[y * stride + x * 4];
                    p[0] = p[1] = p[2] = value;
                    p[3] = 255;
                }
            });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Data, bmData);
            BitmapFormat.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            Bitmap current = BitmapFormat.Bitmap(bmData);
            Bitmap Src = (Bitmap)current.Clone();
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Src, bmSrc);
            Src.Dispose();
            current.Dispose();
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            var Src = (Bitmap)Data.Clone();
            Apply(Data, Src);
            Src.Dispose();
        }
        #endregion

        #region Private methods
        /// <summary>
        /// Apply Gaussian blur.
        /// </summary>
        /// <param name="src">Source</param>
        /// <param name="dst">Destination</param>
        /// <param name="radius">Radius</param>
        private void GaussianBlur(float[,] src, float[,] dst, int radius)
        {
            int size = 2 * radius + 1;
            float[] kernel = new float[size];
            float sigma = radius / 2f;
            float sum = 0;

            for (int i = 0; i < size; i++)
            {
                int x = i - radius;
                kernel[i] = (float)Math.Exp(-(x * x) / (2 * sigma * sigma));
                sum += kernel[i];
            }
            for (int i = 0; i < size; i++)
                kernel[i] /= sum;

            int width = src.GetLength(1);
            int height = src.GetLength(0);

            float[,] temp = new float[height, width];

            // Horizontal blur
            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    float val = 0;
                    for (int k = -radius; k <= radius; k++)
                    {
                        int xx = MathF.Range(x + k, 0, width - 1);
                        val += src[y, xx] * kernel[k + radius];
                    }
                    temp[y, x] = val;
                }
            });

            // Vertical blur
            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    float val = 0;
                    for (int k = -radius; k <= radius; k++)
                    {
                        int yy = MathF.Range(y + k, 0, height - 1);
                        val += temp[yy, x] * kernel[k + radius];
                    }
                    dst[y, x] = val;
                }
            });
        }

        /// <summary>
        /// Apply Sobel filter.
        /// </summary>
        /// <param name="blurred">Blurred source</param>
        /// <param name="gx">Gx</param>
        /// <param name="gy">Gy</param>
        /// <param name="magnitude">Magnitude</param>
        /// <param name="angle">Angle</param>
        private void ApplySobel(float[,] blurred, float[,] gx, float[,] gy, float[,] magnitude, float[,] angle)
        {
            int width = blurred.GetLength(1);
            int height = blurred.GetLength(0);

            int[,] sobelX = new int[,] 
            {
                { -1, 0, 1 },
                { -2, 0, 2 },
                { -1, 0, 1 }
            };

            int[,] sobelY = new int[,] 
            {
                { -1, -2, -1 },
                {  0,  0,  0 },
                {  1,  2,  1 }
            };

            Parallel.For(1, height - 1, y =>
            {
                for (int x = 1; x < width - 1; x++)
                {
                    float sumX = 0, sumY = 0;
                    for (int dy = -1; dy <= 1; dy++)
                        for (int dx = -1; dx <= 1; dx++)
                        {
                            float val = blurred[y + dy, x + dx];
                            sumX += val * sobelX[dy + 1, dx + 1];
                            sumY += val * sobelY[dy + 1, dx + 1];
                        }

                    gx[y, x] = sumX;
                    gy[y, x] = sumY;
                    magnitude[y, x] = MathF.Sqrt(sumX * sumX + sumY * sumY);
                    angle[y, x] = (float)Math.Atan2(sumY, sumX);
                }
            });
        }

        /// <summary>
        /// Non-max suppression.
        /// </summary>
        /// <param name="mag">Magnitude</param>
        /// <param name="dir">Dir</param>
        /// <param name="output">Output</param>
        private void NonMaximumSuppression(float[,] mag, float[,] dir, float[,] output)
        {
            int width = mag.GetLength(1);
            int height = mag.GetLength(0);

            Parallel.For(1, height - 1, y =>
            {
                for (int x = 1; x < width - 1; x++)
                {
                    float angle = dir[y, x] * (180f / MathF.Pi);
                    if (angle < 0) angle += 180f;

                    float current = mag[y, x];
                    float q = 0, r = 0;

                    if ((angle >= 0 && angle < 22.5) || (angle >= 157.5 && angle <= 180))
                    {
                        q = mag[y, x + 1];
                        r = mag[y, x - 1];
                    }
                    else if (angle >= 22.5 && angle < 67.5)
                    {
                        q = mag[y + 1, x - 1];
                        r = mag[y - 1, x + 1];
                    }
                    else if (angle >= 67.5 && angle < 112.5)
                    {
                        q = mag[y + 1, x];
                        r = mag[y - 1, x];
                    }
                    else if (angle >= 112.5 && angle < 157.5)
                    {
                        q = mag[y - 1, x - 1];
                        r = mag[y + 1, x + 1];
                    }

                    output[y, x] = (current >= q && current >= r) ? current : 0;
                }
            });
        }

        /// <summary>
        /// Hysteresis.
        /// </summary>
        /// <param name="nms">NMS</param>
        /// <param name="edges">Edges</param>
        /// <param name="low">Low</param>
        /// <param name="high">High</param>
        private void Hysteresis(float[,] nms, bool[,] edges, float low, float high)
        {
            int width = nms.GetLength(1);
            int height = nms.GetLength(0);

            bool[,] visited = new bool[height, width];

            for (int y = 1; y < height - 1; y++)
            {
                for (int x = 1; x < width - 1; x++)
                {
                    if (!visited[y, x] && nms[y, x] >= high)
                    {
                        FollowEdge(x, y, nms, edges, visited, low);
                    }
                }
            }
        }

        /// <summary>
        /// Follow edge.
        /// </summary>
        /// <param name="x">X</param>
        /// <param name="y">Y</param>
        /// <param name="nms">NMS</param>
        /// <param name="edges">Edges</param>
        /// <param name="visited">Visited</param>
        /// <param name="low">Low</param>
        private void FollowEdge(int x, int y, float[,] nms, bool[,] edges, bool[,] visited, float low)
        {
            if (x < 1 || y < 1 || x >= nms.GetLength(1) - 1 || y >= nms.GetLength(0) - 1)
                return;

            Stack<(int x, int y)> stack = new Stack<(int x, int y)>();
            stack.Push((x, y));

            while (stack.Count > 0)
            {
                (int cx, int cy) = stack.Pop();

                if (visited[cy, cx])
                    continue;

                visited[cy, cx] = true;
                edges[cy, cx] = true;

                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++)
                    {
                        int nx = cx + dx;
                        int ny = cy + dy;

                        if (!visited[ny, nx] && nms[ny, nx] >= low)
                            stack.Push((nx, ny));
                    }
            }
        }
        #endregion
    }

}
