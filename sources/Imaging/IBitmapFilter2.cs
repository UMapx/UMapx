using System.Drawing;
//using System.Drawing.Imaging;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the interface for two images bitmap filter.
    /// </summary>
    public interface IBitmapFilter2
    {
        #region Filter components
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        void Apply(BitmapData bmData, BitmapData bmSrc);
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        void Apply(Bitmap Data, Bitmap Src);
        #endregion
    }
}
