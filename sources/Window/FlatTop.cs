using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the "Flat-Top" window function.
    /// </summary>
    [Serializable]
    public class FlatTop : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the "Flat-Top" window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public FlatTop(int frameSize)
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
            return 1 - 1.93 * Cosine.cosinefunc(2 * x, frameSize) + 1.29 * Cosine.cosinefunc(4 * x, frameSize) - 0.388 * Cosine.cosinefunc(6 * x, frameSize) + 0.028 * Cosine.cosinefunc(8 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            // window function on a discrete time:
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
