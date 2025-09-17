using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the transform base class.
    /// </summary>
    public class TransformBase
    {
        #region Properties
        /// <summary>
        /// Gets or sets the processing direction.
        /// </summary>
        public virtual Direction Direction { get; set; }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public virtual bool Normalized { get; set; }
        #endregion
    }
}
