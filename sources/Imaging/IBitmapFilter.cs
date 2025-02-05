//using System.Drawing;
//using System.Drawing.Imaging;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the interface for bitmap filter.
    /// </summary>
    public interface IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        void Apply(BitmapData bmData);
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        void Apply(Bitmap Data);
        #endregion
    }
}
