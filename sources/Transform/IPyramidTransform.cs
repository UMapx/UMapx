using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the pyramid transform interface.
    /// </summary>
    public interface IPyramidTransform
    {
        #region Interface
        /// <summary>
        /// Forward pyramid transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Pyramid</returns>
        float[][,] Forward(float[,] data);
        /// <summary>
        /// Backward pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        float[,] Backward(float[][,] pyramid);
        /// <summary>
        /// Forward pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        float[][] Forward(float[] data);
        /// <summary>
        /// Backward pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        float[] Backward(float[][] pyramid);
        /// <summary>
        /// Forward pyramid transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Pyramid</returns>
        ComplexF[][,] Forward(ComplexF[,] data);
        /// <summary>
        /// Backward pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        ComplexF[,] Backward(ComplexF[][,] pyramid);
        /// <summary>
        /// Forward pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        ComplexF[][] Forward(ComplexF[] data);
        /// <summary>
        /// Backward pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        ComplexF[] Backward(ComplexF[][] pyramid);
        #endregion
    }
}
