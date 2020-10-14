using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Lanczos window function.
    /// </summary>
    [Serializable]
    public class Lanzcos : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Lanczos window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Lanzcos(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            // Lanczos function:
            return Special.Sinc(2 * x / (frameSize - 1) - 1, Math.PI);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
