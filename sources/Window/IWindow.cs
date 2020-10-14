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
        /// <returns>Double precision floating point number</returns>
        double Function(double x);
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        double[] GetWindow();
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        double[] GetWindow(int frameSize);
        /// <summary>
        /// Returns an array of window function values.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        double[] Function(double[] x, int frameSize);
        /// <summary>
        /// Returns an array of window function values.
        /// </summary>
        /// <param name="x">Array</param>
        /// <returns>Array</returns>
        double[] Function(double[] x);
        #endregion
    }
}
