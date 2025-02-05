using System;
using System.Drawing;
using UMapx.Core;

//using System.Drawing.Imaging;
using SkiaDrawing;


namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the additive noise filter.
    /// <remarks>
    /// Filter usage example:
    /// https://en.wikipedia.org/wiki/Gaussian_noise
    /// </remarks>
    /// </summary>
    [Serializable]
    public class AdditiveNoise : IBitmapFilter
    {
        #region Private data
        private int amount = 10;
        private Random generator = new Random();
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the additive noise filter.
        /// </summary>
        public AdditiveNoise() { }
        /// <summary>
        /// Initializes the additive noise filter.
        /// </summary>
        /// <param name="amount">Amount [0, 100]</param>
        public AdditiveNoise(int amount)
        {
            Amount = amount;
        }
        /// <summary>
        /// Gets or sets the amout value [0, 100].
        /// </summary>
        public int Amount
        {
            get
            {
                return this.amount;
            }
            set
            {
                this.amount = Math.Max(0, Math.Min(100, value));
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        //public unsafe void Apply(BitmapData bmData)
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            int stride = bmData.Stride;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[2] = Maths.Byte(p[2] + generator.Next(-amount, amount));
                    p[1] = Maths.Byte(p[1] + generator.Next(-amount, amount));
                    p[0] = Maths.Byte(p[0] + generator.Next(-amount, amount));
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
    }
}
