using System;
using System.Drawing.Imaging;
using System.Drawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines a normal bump filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Normal_mapping
    /// </remarks>
    /// </summary>
    [Serializable]
    public class NormalBump : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private readonly RGBFilter rgbFilter = new RGBFilter(0, 0, 255);
        private Convolution convolutionFilter;
        private int scale = 1;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the normal bump filter.
        /// </summary>
        /// <param name="scale">Scale [-5, 5]</param>
        public NormalBump(int scale = 0)
        {
            Scale = scale;
        }
        /// <summary>
        /// Gets or sets the scale value [-5, 5].
        /// </summary>
        public int Scale
        { 
            get
            {
                return scale;
            }
            set
            {
                this.scale = value;
                this.convolutionFilter = new Convolution(NormalBumpOperator(this.scale), 128, false);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            convolutionFilter.Apply(bmData, bmSrc);
            rgbFilter.Apply(bmData);
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

        #region Public static
        /// <summary>
        /// Returns the normal bump operator.
        /// </summary>
        /// <param name="scale">Scale</param>
        /// <returns>Matrix</returns>
        public static float[,] NormalBumpOperator(int scale)
        {
            var T0 = +1 + scale;
            var T1 = -1 - scale * 2;

            return new float[,] { { T0, -1, T0 }, { -1, 0, -1 }, { 1, T1, 1 } };
        }
        #endregion
    }
}
