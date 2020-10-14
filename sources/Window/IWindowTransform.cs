namespace UMapx.Window
{
    /// <summary>
    /// Defines the general window transform interface.
    /// </summary>
    public interface IWindowTransform
    {
        #region Interface
        /// <summary>
        /// Gets or sets the window function.
        /// </summary>
        IWindow Window { get; set; }
        #endregion
    }
}
