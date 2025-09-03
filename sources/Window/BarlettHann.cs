using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Barlett-Hann window function.
    /// </summary>
    [Serializable]
    public class BarlettHann : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Barlett-Hann window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public BarlettHann(int frameSize)
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
            float a = Math.Abs(Math.Abs(x / (frameSize - 1)) - 0.5f);
            return 0.62f - 0.48f * a - 0.38f * Cosine.Cosinefunc(2 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            float t = frameSize - 1;
            float[] x = MatrixF.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
