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
        /// <returns>Value</returns>
        public override float Function(float x, int frameSize)
        {
            // Lanczos function:
            return Special.Sinc(2 * x / (frameSize - 1) - 1, Maths.Pi);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            float t = (frameSize - 1);
            float[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
