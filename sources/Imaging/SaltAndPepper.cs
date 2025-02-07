using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the salt and pepper noise filter.
    /// <remarks>
    /// Filter usage example:
    /// https://en.wikipedia.org/wiki/Salt-and-pepper_noise
    /// </remarks>
    /// </summary>
    [Serializable]
    public class SaltAndPepper : IBitmapFilter
    {
        #region Private data
        private double amount = 10;
        private Random generator = new Random();
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the salt and pepper noise filter.
        /// </summary>
        public SaltAndPepper() { }
        /// <summary>
        /// Initializes the salt and pepper noise filter.
        /// </summary>
        /// <param name="amount">Amount [0, 100].</param>
        public SaltAndPepper(double amount)
        {
            Amount = amount;
        }
        /// <summary>
        /// Gets or sets the amout value [0, 100].
        /// </summary>
        public double Amount
        {
            get
            {
                return amount;
            }
            set
            {
                amount = Math.Max(0, Math.Min(100, value));
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            int noisyPixels = (int)((width * height * amount) / 100);
            byte[] values = new byte[2] { 0, 255 };
            int stride = bmData.Stride;
            int i, colorPlane;

            for (i = 0; i < noisyPixels; i++)
            {
                x = generator.Next(width);
                y = generator.Next(height);
                colorPlane = generator.Next(3);
                p[y * stride + x * 4 + colorPlane] = values[generator.Next(2)];
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
    }
}
