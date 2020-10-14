using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Gabor window function.
    /// </summary>
    [Serializable]
    public class Gabor : WindowBase
    {
        #region Private data
        private double sigma = 1;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Gabor window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="sigma">Scale parameter</param>
        public Gabor(int frameSize, double sigma = 1)
        {
            this.FrameSize = frameSize;
            this.Sigma = sigma;
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
            // Gabor window function
            double y = x / frameSize;
            double z = Math.Pow(2 * Math.PI * y / sigma, 2);
            return Math.Exp(-z);
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
    }
}
