using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Kaiser window function.
    /// </summary>
    [Serializable]
    public class Kaiser : WindowBase
    {
        #region Private data
        private double a = 3;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Kaiser window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="a">Form parameter</param>
        public Kaiser(int frameSize, double a = 3)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Gets or sets the value of the form parameter.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
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
            // Kaiser window:
            double u = 2 * x / (frameSize - 1);
            double r = 1 - u * u;
            double v = r >= 0 ? Math.Sqrt(1 - u * u) : 0;
            double z = Math.PI * this.a;
            double q = Special.I(z * v, 0);
            return q / Special.I(z, 0);
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
