using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Hamming window function.
    /// </summary>
    [Serializable]
    public class Hamming : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Hamming window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Hamming(int frameSize)
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
            return 0.53836f - 0.46164f * Cosine.Cosinefunc(2 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            float t = frameSize - 1;
            float[] x = Matrix.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
