using System;
using System.Collections.Generic;
using System.Text;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the cosine window function.
    /// </summary>
    [Serializable]
    public class Cosine : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the cosine window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Cosine(int frameSize)
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
            return Math.Cos(Math.PI * x / (frameSize - 1));
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion

        #region Static components
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Factor</returns>
        internal static double cosinefunc(double x, int frameSize)
        {
            return Math.Cos(Math.PI * x / (frameSize - 1));
        }
        #endregion
    }
}
