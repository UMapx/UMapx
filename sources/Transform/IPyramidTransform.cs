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
        double[][,] Forward(double[,] data);
        /// <summary>
        /// Backward pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        double[,] Backward(double[][,] pyramid);
        /// <summary>
        /// Forward pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        double[][] Forward(double[] data);
        /// <summary>
        /// Backward pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        double[] Backward(double[][] pyramid);
        /// <summary>
        /// Forward pyramid transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Pyramid</returns>
        Complex[][,] Forward(Complex[,] data);
        /// <summary>
        /// Backward pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        Complex[,] Backward(Complex[][,] pyramid);
        /// <summary>
        /// Forward pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        Complex[][] Forward(Complex[] data);
        /// <summary>
        /// Backward pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        Complex[] Backward(Complex[][] pyramid);
        #endregion
    }
}
