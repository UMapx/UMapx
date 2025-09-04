using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Welch window function.
    /// </summary>
    [Serializable]
    public class Welch : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Welch window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Welch(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Value</returns>
        public override float Function(float x, int frameSize)
        {
            // Welch function:
            float t = (frameSize - 1) / 2.0f;
            float a = (x - t) / t;
            return 1 - a * a;
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
