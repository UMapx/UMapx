using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Colorspace;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the vibrance correction filter (RGB-based).
    /// </summary>
    [Serializable]
    public class VibranceCorrection : IBitmapFilter
    {
        #region Private data
        private float vibrance;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the vibrance correction filter.
        /// </summary>
        /// <param name="vibrance">Vibrance [-100, 100]</param>
        public VibranceCorrection(float vibrance)
        {
            Vibrance = vibrance;
        }

        /// <summary>
        /// Initializes the vibrance correction filter.
        /// </summary>
        public VibranceCorrection()
        {
            Vibrance = 20;
        }

        /// <summary>
        /// Gets or sets the vibrance value [-100, 100].
        /// </summary>
        public float Vibrance
        {
            get => vibrance;
            set => vibrance = Math.Max(-100, Math.Min(100, value));
        }

        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            float v = this.vibrance / 100.0f;

            Parallel.For(0, height, y =>
            {
                int ystride = y * stride;

                for (int x = 0; x < width; x++)
                {
                    int k = ystride + x * 4;

                    byte r = p[k + 2];
                    byte g = p[k + 1];
                    byte b = p[k + 0];
                    
                    RGB rgb = RGB.Vibrance(r, g, b, v);

                    p[k + 2] = rgb.Red;
                    p[k + 1] = rgb.Green;
                    p[k + 0] = rgb.Blue;
                }
            });
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
