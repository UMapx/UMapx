using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Zak transform interface.
    /// </summary>
    public interface IZakTransform
    {
        #region Interface
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Matrix</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        ComplexF[,] Forward(ComplexF[] input);
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Matrix</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        ComplexF[,] Forward(float[] input);
        /// <summary>
        /// Backward Zak transform.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <returns>Array</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        ComplexF[] Backward(ComplexF[,] matrix);
        #endregion
    }
}
