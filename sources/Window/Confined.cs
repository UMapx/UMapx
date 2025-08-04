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
        private float sigma = 1;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the closed Gaussian window.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="sigma">Standard deviation (0.14 * N)</param>
        public Confined(int frameSize, float sigma = 1)
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
            this.Sigma = 0.14f * frameSize;
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Gets or sets the standard deviation (>0).
        /// </summary>
        public float Sigma
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
        /// <returns>Value</returns>
        public override float Function(float x, int frameSize)
        {
            float a = G(-0.5f) * (G(x + frameSize) + G(x - frameSize));
            float b = G(-0.5f + frameSize) + G(-0.5f - frameSize);
            return G(x) - a / b;
        }
        /// <summary>
        /// Функция G(x).
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Value</returns>
        private float G(float x)
        {
            float a = (frameSize - 1) / 2;
            float t = (x - a) / (2 * sigma);
            return (float)Math.Exp(-t * t);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            // window function on a discrete time:
            float t = frameSize - 1;
            float[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
