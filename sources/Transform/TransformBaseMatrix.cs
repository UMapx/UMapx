namespace UMapx.Transform
{
    /// <summary>
    /// Defines an adapter class for a matrix transform.
    /// </summary>
    public abstract class TransformBaseMatrix : TransformBase
    {
        #region Special methods
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
