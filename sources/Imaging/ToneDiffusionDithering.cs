using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the tone diffusion dithering filter.
    /// <remarks>
    /// Filter usage example:
    /// https://en.wikipedia.org/wiki/Dither
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ToneDiffusionDithering : IBitmapFilter
    {
        #region Private data
        private double[,] matrix;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the tone diffusion dithering filter.
        /// </summary>
        public ToneDiffusionDithering()
        {
            this.matrix = new double[4, 4] {
                                        {  15, 143,  47, 175 },
                                        { 207,  79, 239, 111 },
                                        {  63, 191,  31, 159 },
                                        { 255, 127, 223,  95 }};
        }
        /// <summary>
        /// Initializes the tone diffusion dithering filter.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        public ToneDiffusionDithering(double[,] matrix)
        {
            Matrix = matrix;
        }
        /// <summary>
        /// Gets or sets the tone diffusion dithering matrix. 
        /// </summary>
        public double[,] Matrix
        {
            get
            {
                return this.matrix;
            }
            set
            {
                this.matrix = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int rows = matrix.GetLength(0), cols = matrix.GetLength(1);
            int y, x, width = bmData.Width, height = bmData.Height;
            byte n;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    n = (byte)(matrix[(y % rows), (x % cols)]);
                    p[2] = (byte)((p[2] <= n) ? 0 : 255);
                    p[1] = (byte)((p[1] <= n) ? 0 : 255);
                    p[0] = (byte)((p[0] <= n) ? 0 : 255);
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            Apply(bmData);
            BitmapFormat.Unlock(Data, bmData);
        }
        #endregion

        #region Static methods
        private static System.Random rand = new System.Random();
        /// <summary>
        /// Initializes the order dithering filter.
        /// <remarks>
        /// More information can be found on the website:
        /// http://en.wikipedia.org/wiki/Ordered_dithering
        /// Filter usage example:
        /// https://en.wikipedia.org/wiki/Dither
        /// </remarks>
        /// </summary>
        /// <param name="radius">Radius [0, 255]</param>
        /// <returns>Tone diffusion dithering filter</returns>
        public static ToneDiffusionDithering Order(int radius)
        {
            byte c = (byte)(256 / radius / radius + 1), d = 0;
            double[,] table = new double[radius, radius];
            int i, j;

            for (i = 0; i < radius; i++)
            {
                for (j = 0; j < radius; j++, d += c)
                {
                    table[i, j] = d;
                }
            }

            return new ToneDiffusionDithering(table);
        }
        /// <summary>
        /// Initializes the random dithering filter.
        /// </summary>
        /// <param name="radius">Radius [0, 255]</param>
        /// <returns>Tone diffusion dithering filter</returns>
        public static ToneDiffusionDithering Random(int radius)
        {
            double[,] table = new double[radius, radius];
            int i, j;

            for (i = 0; i < radius; i++)
            {
                for (j = 0; j < radius; j++)
                {
                    table[i, j] = (byte)rand.Next(0, 255);
                }
            }

            return new ToneDiffusionDithering(table);
        }
        /// <summary>
        /// Initializes the classic dithering filter.
        /// </summary>
        /// <returns>Tone diffusion dithering filter</returns>
        public static ToneDiffusionDithering Basic()
        {
            return new ToneDiffusionDithering(new double[4, 4] {
                                        {  15, 143,  47, 175 },
                                        { 207,  79, 239, 111 },
                                        {  63, 191,  31, 159 },
                                        { 255, 127, 223,  95 }});
        }
        /// <summary>
        /// Initializes the Bayer dithering filter.
        /// </summary>
        /// <returns>Tone diffusion dithering filter</returns>
        public static ToneDiffusionDithering Bayer()
        {
            return new ToneDiffusionDithering(new double[,] {
                                      {   0, 192,  48, 240 },
                                      { 128,  64, 176, 112 },
                                      {  32, 224,  16, 208 },
                                      { 160,  96, 144,  80 } });
        }
        #endregion
    }
}
