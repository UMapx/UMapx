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
        private float a = 3;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Kaiser window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="a">Form parameter</param>
        public Kaiser(int frameSize, float a = 3)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Gets or sets the value of the form parameter.
        /// </summary>
        public float A
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
        /// <param name="x">Value</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Value</returns>
        public override float Function(float x, int frameSize)
        {
            // Kaiser window:
            float u = 2 * x / (frameSize - 1);
            float r = 1 - u * u;
            float v = r >= 0 ? (float)Math.Sqrt(1 - u * u) : 0;
            float z = Maths.Pi * this.a;
            float q = Special.I(z * v, 0);
            return q / Special.I(z, 0);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            float t = (frameSize - 1) / 2.0f;
            float[] x = Matrix.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
