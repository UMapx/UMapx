using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Welch window function.
    /// </summary>
    [Serializable]
    public class Welch : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Welch window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Welch(int frameSize)
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
            // Welch function:
            double t = (frameSize - 1) / 2.0;
            double a = (x - t) / t;
            return 1 - a * a;
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
