using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the transform interface.
    /// </summary>
    public interface ITransform
    {
        #region Interface
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Array</returns>
        float[] Forward(float[] data);
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Array</returns>
        float[,] Forward(float[,] data);
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Array</returns>
        float[] Backward(float[] data);
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Matrix</returns>
        float[,] Backward(float[,] data);
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Array</returns>
        Complex[] Forward(Complex[] data);
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Array</returns>
        Complex[,] Forward(Complex[,] data);
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Array</returns>
        Complex[] Backward(Complex[] data);
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Matrix</returns>
        Complex[,] Backward(Complex[,] data);
        #endregion
    }
}
