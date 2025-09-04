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
        private float a = 3;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Tukey window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="a">Form parameter [0, 1]</param>
        public Tukey(int frameSize, float a = 1)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Gets or sets the value of the form parameter [0, 1].
        /// </summary>
        public float A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = Maths.Float(value);
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
            // Tukey window:
            float n = frameSize - 1;
            float d = n * (1 - a / 2.0f);
            float b = n / 2.0f;
            float c = a * b;

            // Creating:
            if (x >= 0 && x < c)
            {
                return 0.5f * (1 + (float)Math.Cos(Math.PI * (x / c - 1)));
            }
            else if (x >= c && x <= d)
            {
                return 1.0f;
            }
            else if (x > d && x <= n)
            {
                return 0.5f * (1 + (float)Math.Cos(Math.PI * (x / c - 2.0 / a + 1)));
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
            float[] x = Matrix.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
