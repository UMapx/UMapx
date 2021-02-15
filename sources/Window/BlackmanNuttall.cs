using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Blackman-Nuttall window function.
    /// </summary>
    [Serializable]
    public class BlackmanNuttall : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Blackman-Nuttall window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public BlackmanNuttall(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>float precision floating point number</returns>
        public override float Function(float x, int frameSize)
        {
            return 0.3635819f - 0.4891775f * Cosine.cosinefunc(2 * x, frameSize) + 0.1365995f * Cosine.cosinefunc(4 * x, frameSize) - 0.0106411f * Cosine.cosinefunc(6 * x, frameSize);
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
