using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the sine window function.
    /// </summary>
    [Serializable]
    public class Sine : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the sine window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Sine(int frameSize)
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
            return Math.Sin(Math.PI * x / (frameSize - 1));
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

        #region Static components
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        internal static double sinefunc(double x, int frameSize)
        {
            return Math.Sin(Math.PI * x / (frameSize - 1));
        }
        #endregion
    }
}
