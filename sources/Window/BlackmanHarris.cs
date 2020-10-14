using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Blackman-Harris window function.
    /// </summary>
    [Serializable]
    public class BlackmanHarris : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Blackman-Harris window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public BlackmanHarris(int frameSize)
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
            return 0.35875 - 0.48829 * Cosine.cosinefunc(2 * x, frameSize) + 0.14128 * Cosine.cosinefunc(4 * x, frameSize) - 0.01168 * Cosine.cosinefunc(6 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
