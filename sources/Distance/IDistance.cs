using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines the distance interface.
    /// </summary>
    internal interface IDistance
    {
        #region Interface
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        float Compute(float[] p, float[] q);
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        Complex32 Compute(Complex32[] p, Complex32[] q);
        /// <summary>
        /// Returns distance values. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        float[] Compute(float[,] p, float[,] q);
        /// <summary>
        /// Returns distance values. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        Complex32[] Compute(Complex32[,] p, Complex32[,] q);
        #endregion
    }
}
