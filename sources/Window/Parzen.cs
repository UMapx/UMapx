using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Parzen window function.
    /// </summary>
    [Serializable]
    public class Parzen : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Parzen window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Parzen(int frameSize)
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
            // coefficients:
            double y = Math.Abs(x);
            double c = frameSize / 2.0;
            double a = x / c;
            double b = y / c;

            // props:
            if ((y >= 0) &&
                (y <= frameSize / 4))
            {
                return 1.0 - 6.0 * a * a * (1.0 - b);
            }
            else if (y >= frameSize / 4 &&
                y <= frameSize / 2)
            {
                return 2 * Math.Pow(1 - b, 3);
            }
            return 0.0;
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
