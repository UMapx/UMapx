using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the interface of canvas.
    /// </summary>
    public interface ICanvas
    {
        #region Interface
        /// <summary>
        /// Gets or sets the width of the canvas.
        /// </summary>
        int Width { get; set; }
        /// <summary>
        /// Gets or sets the height of the canvas.
        /// </summary>
        int Height { get; set; }
        /// <summary>
        /// Creates canvas.
        /// </summary>
        /// <returns>Bitmap</returns>
        Bitmap Create();
        #endregion
    }
}
