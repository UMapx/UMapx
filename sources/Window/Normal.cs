using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines a generalized window normal function.
    /// </summary>
    [Serializable]
    public class Normal : WindowBase
    {
        #region Private data
        private double sigma = 1;
        private double p = 2;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes a generalized window normal function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="sigma">Standard deviation (>0)</param>
        /// <param name="pow">Power<remarks>For p = 2 - Gaussian window</remarks></param>
        public Normal(int frameSize, double sigma = 1, double pow = 2)
        {
            this.Sigma = sigma;
            this.FrameSize = frameSize;
            this.p = pow;
        }
        /// <summary>
        /// Gets or sets the standard deviation (>0).
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Power.
        /// </summary>
        public double Pow
        {
            get
            {
                return this.p;
            }
            set
            {
                this.p = value;
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
            double a = (frameSize - 1) / 2;
            double t = (x - a) / (sigma * a);
            return Math.Exp(-Math.Pow(t, p));
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
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
