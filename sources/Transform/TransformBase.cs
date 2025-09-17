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

        #region Normalization methods
        /// <summary>
        /// Returns the forward normalization factor.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Value</returns>
        protected virtual float ForwardNormalizationFactor(int n)
        {
            return 1.0f;
        }
        /// <summary>
        /// Returns the backward normalization factor.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Value</returns>
        protected virtual float BackwardNormalizationFactor(int n)
        {
            return 1.0f;
        }
        /// <summary>
        /// Returns the forward matrix size.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Value</returns>
        protected virtual int ForwardMatrixSize(int n)
        {
            return n;
        }
        /// <summary>
        /// Returns the backward matrix size.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Value</returns>
        protected virtual int BackwardMatrixSize(int n)
        {
            return n;
        }
        #endregion
    }
}
