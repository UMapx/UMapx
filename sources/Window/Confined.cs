using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the closed Gaussian window.
    /// </summary>
    [Serializable]
    public class Confined : WindowBase
    {
        #region Private data
        private double sigma = 1;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the closed Gaussian window.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="sigma">Standard deviation (0.14 * N)</param>
        public Confined(int frameSize, double sigma = 1)
        {
            this.Sigma = sigma;
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Initializes a Gaussian window function closed.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Confined(int frameSize)
        {
            this.Sigma = 0.14 * frameSize;
            this.FrameSize = frameSize;
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
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            // Вычисление функции:
            double a = G(-0.5) * (G(x + frameSize) + G(x - frameSize));
            double b = G(-0.5 + frameSize) + G(-0.5 - frameSize);
            return G(x) - a / b;
        }
        /// <summary>
        /// Функция G(x).
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        private double G(double x)
        {
            double a = (frameSize - 1) / 2;
            double t = (x - a) / (2 * sigma);
            return Math.Exp(-t * t);
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
