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
        private float sigma = 1;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Gabor window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="sigma">Scale parameter</param>
        public Gabor(int frameSize, float sigma = 1)
        {
            this.FrameSize = frameSize;
            this.Sigma = sigma;
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
        /// <param name="x">Value</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Value</returns>
        public override float Function(float x, int frameSize)
        {
            // Gabor window function
            float y = x / frameSize;
            float z = (float)Math.Pow(2 * Math.PI * y / sigma, 2);
            return (float)Math.Exp(-z);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            float t = (frameSize - 1) / 2.0f;
            float[] x = MatrixF.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion

        #region Static methods
        /// <summary>
        /// Returns Gabor window function defined without sigma. Scaled function version.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Gabor window function</returns>
        public static Gabor Scaled(int frameSize)
        {
            var s = 1.0f / MathF.Sqrt(frameSize / 4.0f);
            return new Gabor(frameSize, s);
        }
        #endregion
    }
}
