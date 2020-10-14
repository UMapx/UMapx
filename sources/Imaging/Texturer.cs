using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the texturing filter.
    /// </summary>
    [Serializable]
    public class Texturer : IBitmapFilter
    {
        #region Private data
        private double[,] texture = null;
        private double depth = 0.5;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the texturing filter.
        /// </summary>
        /// <param name="texture">Matrix</param>
        public Texturer(double[,] texture)
        {
            Texture = texture;
        }
        /// <summary>
        /// Initializes the texturing filter.
        /// </summary>
        /// <param name="texture">Matrix</param>
        /// <param name="depth">Depth [0, 1]</param>
        public Texturer(double[,] texture, double depth = 1.0)
        {
            Texture = texture; Depth = depth;
        }
        /// <summary>
        /// Gets or sets the texture matrix.
        /// </summary>
        public double[,] Texture
        {
            get { return texture; }
            set { texture = value; }
        }
        /// <summary>
        /// Gets or sets the depth value [0, 1].
        /// </summary>
        public double Depth
        {
            get { return this.depth; }
            set { this.depth = Maths.Double(value); }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            int widthToProcess = Math.Min(width, texture.GetLength(1));
            int heightToProcess = Math.Min(height, texture.GetLength(0));
            double z = 1.0 - this.depth;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, heightToProcess, y =>
            {
                int x, ystride, k;
                byte red, green, blue;
                double t;

                ystride = y * stride;

                for (x = 0; x < widthToProcess; x++)
                {
                    k = ystride + x * 4;
                    red = p[k + 2]; green = p[k + 1]; blue = p[k];
                    t = texture[y, x];

                    p[k + 2] = Maths.Byte((z * red) + (this.depth * red) * t);
                    p[k + 1] = Maths.Byte((z * green) + (this.depth * green) * t);
                    p[k] = Maths.Byte((z * blue) + (this.depth * blue) * t);
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Static methods
        private static Random rand = new Random();
        private static int r;

        /// <summary>
        /// Implements the construction of a wood texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="rings">Rings</param>
        /// <returns>Matrix</returns>
        public static Texturer Wood(int m, int l, double rings = 12)
        {
            PerlinNoise noise = new PerlinNoise(8, 0.5, 1.0 / 32, 0.05); r = rand.Next(5000);
            double[,] texture = new double[m, l];
            int w2 = l / 2, h2 = m / 2;

            Parallel.For(0, m, y =>
            {
                double xv, yv;
                int x;

                for (x = 0; x < l; x++)
                {
                    xv = (double)(x - w2) / l;
                    yv = (double)(y - h2) / m;

                    texture[y, x] = Math.Min(1.0f, (float)Math.Abs(Math.Sin((Math.Sqrt(xv * xv + yv * yv) + noise.Function2D(x + r, y + r)) * Math.PI * 2 * rings)));
                }
            }
            );

            return new Texturer(texture);
        }
        /// <summary>
        /// Implements the construction of a textile texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Texturer Textile(int m, int l)
        {
            PerlinNoise noise = new PerlinNoise(3, 0.65, 1.0 / 8, 1.0); r = rand.Next(5000);
            double[,] texture = new double[m, l];

            Parallel.For(0, m, y =>
            {
                int x;
                for (x = 0; x < l; x++)
                {
                    texture[y, x] = Math.Max(0.0f, Math.Min(1.0f, (
                                (double)Math.Sin(x + noise.Function2D(x + r, y + r)) +
                                (double)Math.Sin(y + noise.Function2D(x + r, y + r))) * 0.25f + 0.5f));
                }
            }
            );

            return new Texturer(texture);
        }
        /// <summary>
        /// Implements the construction of a marble texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="yPeriod">Y-period</param>
        /// <param name="xPeriod">X-period</param>
        /// <returns>Matrix</returns>
        public static Texturer Marble(int m, int l, double yPeriod = 10.0, double xPeriod = 5.0)
        {
            PerlinNoise noise = new PerlinNoise(2, 0.65, 1.0 / 32, 1.0); r = rand.Next(5000);
            double[,] texture = new double[m, l];
            double xFact = xPeriod / l;
            double yFact = yPeriod / m;

            Parallel.For(0, m, y =>
            {
                int x;

                for (x = 0; x < l; x++)
                {
                    texture[y, x] = Math.Min(1.0f, (float)Math.Abs(Math.Sin((x * xFact + y * yFact + noise.Function2D(x + r, y + r)) * Math.PI)));
                }
            }
            );

            return new Texturer(texture);
        }
        /// <summary>
        /// Implements the construction of a labyrinth texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Texturer Labyrinth(int m, int l)
        {
            PerlinNoise noise = new PerlinNoise(1, 0.65, 1.0 / 16, 1.0); r = rand.Next(5000);
            double[,] texture = new double[m, l];

            Parallel.For(0, m, y =>
            {
                int x;

                for (x = 0; x < l; x++)
                {
                    texture[y, x] = Math.Min(1.0f, (float)Math.Abs(noise.Function2D(x + r, y + r)));

                }
            }
            );

            return new Texturer(texture);
        }
        /// <summary>
        /// Implements the construction of a clouds texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Texturer Clouds(int m, int l)
        {
            PerlinNoise noise = new PerlinNoise(8, 0.5, 1.0 / 32, 1.0); r = rand.Next(5000);
            double[,] texture = new double[m, l];

            Parallel.For(0, m, y =>
            {
                int x;

                for (x = 0; x < l; x++)
                {
                    texture[y, x] = Math.Max(0.0f, Math.Min(1.0f, (double)noise.Function2D(x + r, y + r) * 0.5f + 0.5f));

                }
            }
            );

            return new Texturer(texture);
        }
        #endregion
    }
}
