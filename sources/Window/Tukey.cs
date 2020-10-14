using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the window function of Tukey.
    /// </summary>
    [Serializable]
    public class Tukey : WindowBase
    {
        #region Private data
        private double a = 3;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Tukey window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="a">Form parameter [0, 1]</param>
        public Tukey(int frameSize, double a = 1)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Gets or sets the value of the form parameter [0, 1].
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = Maths.Double(value);
            }
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            // Tukey window:
            double n = frameSize - 1;
            double d = n * (1 - a / 2.0);
            double b = n / 2.0;
            double c = a * b;

            // Creating:
            if (x >= 0 && x < c)
            {
                return 0.5 * (1 + Math.Cos(Math.PI * (x / c - 1)));
            }
            else if (x >= c && x <= d)
            {
                return 1.0;
            }
            else if (x > d && x <= n)
            {
                return 0.5 * (1 + Math.Cos(Math.PI * (x / c - 2.0 / a + 1)));
            }
            return 0;
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
