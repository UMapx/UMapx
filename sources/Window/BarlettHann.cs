using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Barlett-Hann window function.
    /// </summary>
    [Serializable]
    public class BarlettHann : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Barlett-Hann window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public BarlettHann(int frameSize)
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
            // Berlett-Hann function:
            double a = Math.Abs(Math.Abs(x / (frameSize - 1)) - 0.5);
            return 0.62 - 0.48 * a - 0.38 * Cosine.cosinefunc(2 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
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
