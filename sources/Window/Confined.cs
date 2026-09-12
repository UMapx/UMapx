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
        /// <param name="frameSize">Window size.</param>
        /// <param name="sigma">Standard deviation (0.14 * N).</param>
        public Confined(int frameSize, float sigma = 1)
        {
            this.Sigma = sigma;
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Initializes a Gaussian window function closed.
        /// </summary>
        /// <param name="frameSize">Window size.</param>
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
                    throw new ArgumentException("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="frameSize">Window size.</param>
        /// <returns>Value.</returns>
        public override float Function(float x, int frameSize)
        {
            if (frameSize == 1) return 1;
            // G(-1/2) equals G(N-1/2) by symmetry. Divide before
            // multiplying to avoid a 0/0 ratio for narrow Gaussians.
            double scale = 2.0 * sigma;
            double edge = frameSize / (2.0 * scale);
            double far = 3.0 * frameSize / (2.0 * scale);
            double ratio = 1.0 / (1.0 + Math.Exp(edge * edge - far * far));
            return (float)(G(x, frameSize) - ratio * (G(x + frameSize, frameSize) + G(x - frameSize, frameSize)));
        }
        /// <summary>
        /// Evaluates the Gaussian centered on the supplied window.
        /// </summary>
        /// <param name="x">Sample coordinate, including points outside the window.</param>
        /// <param name="frameSize">Number of samples defining the center.</param>
        /// <returns>The Gaussian value with scale twice the positive standard deviation.</returns>
        private double G(float x, int frameSize)
        {
            double a = (frameSize - 1) / 2.0;
            double t = (x - a) / (2.0 * sigma);
            return Math.Exp(-t * t);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array.</returns>
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
