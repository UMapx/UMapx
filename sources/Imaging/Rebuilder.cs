namespace UMapx.Imaging
{
    /// <summary>
    /// Defines an abstract data rebuilding class.
    /// </summary>
    public abstract class Rebuilder
    {
        #region Protected components
        /// <summary>
        /// Use data rebuilding or not.
        /// </summary>
        protected bool rebuild = false;
        /// <summary>
        /// Implements the rebuilding of class data.
        /// </summary>
        protected abstract void Rebuild();
        #endregion
    }
}
