namespace UMapx.Window
{
    /// <summary>
    /// Defines the interface of window functions.
    /// </summary>
    public interface IWindow
    {
        #region Interface
        /// <summary>
        /// Gets or sets the window size.
        /// </summary>
        int FrameSize { get; set; }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Value</returns>
        float Function(float x);
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        float[] GetWindow();
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        float[] GetWindow(int frameSize);
        /// <summary>
        /// Returns an array of window function values.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        float[] Function(float[] x, int frameSize);
        /// <summary>
        /// Returns an array of window function values.
        /// </summary>
        /// <param name="x">Array</param>
        /// <returns>Array</returns>
        float[] Function(float[] x);
        #endregion
    }
}
