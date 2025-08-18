using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Zak transform interface.
    /// </summary>
    public interface IZakTransform
    {
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Matrix</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        Complex32[,] Forward(Complex32[] input);
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Matrix</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        Complex32[,] Forward(float[] input);
        /// <summary>
        /// Backward Zak transform.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <returns>Array</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        Complex32[] Backward(Complex32[,] matrix);
    }
}
