using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the window function of Planck.
    /// </summary>
    [Serializable]
    public class Planck : WindowBase
    {
        #region Private data
        private float a = 0.15f;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Planck window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="a">Form parameter [0, 0.5]</param>
        public Planck(int frameSize, float a = 0.15f)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Gets or sets the value of the form parameter [0, 0.5].
        /// </summary>
        public float A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = MathF.Range(value, 0, 0.5f);
            }
        }
        /// <summary>
        /// Function Z+-(x, a).
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="p">Sign</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Value</returns>
        private float Z(float x, bool p, int frameSize)
        {
            // params:
            float t = p ? 1 : -1;
            float y = 2 * x / (frameSize - 1) - 1;

            // function:
            float u = 1.0f / (1 + t * y);
            float v = 1.0f / (1 - 2 * a + t * y);
            return 2 * a * (u + v);
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Value</returns>
        public override float Function(float x, int frameSize)
        {
            // Planck taper window:
            float n = frameSize - 1;
            float b = a * n;
            float c = (1 - a) * n;

            // Creating:
            if (x >= 0 && x < b)
            {
                return 1.0f / (float)(Math.Exp(Z(x, true, frameSize)) + 1);
            }
            else if (x >= b && x <= c)
            {
                return 1.0f;
            }
            else if (x > c && x <= n)
            {
                return 1.0f / (float)(Math.Exp(Z(x, false, frameSize)) + 1);
            }
            return 0;
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            float t = (frameSize - 1);
            float[] x = MatrixF.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
