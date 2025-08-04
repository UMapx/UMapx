using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the "Flat-Top" window function.
    /// </summary>
    [Serializable]
    public class FlatTop : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the "Flat-Top" window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public FlatTop(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Value</returns>
        public override float Function(float x, int frameSize)
        {
            return 1 - 1.93f * Cosine.cosinefunc(2 * x, frameSize) + 1.29f * Cosine.cosinefunc(4 * x, frameSize) - 0.388f * Cosine.cosinefunc(6 * x, frameSize) + 0.028f * Cosine.cosinefunc(8 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
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
