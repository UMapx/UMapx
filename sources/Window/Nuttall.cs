using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Nuttall window function.
    /// </summary>
    [Serializable]
    public class Nuttall : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Nuttall window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Nuttall(int frameSize)
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
            return 0.355768f - 0.487396f * Cosine.cosinefunc(2 * x, frameSize) + 0.144232f * Cosine.cosinefunc(4 * x, frameSize) - 0.012604f * Cosine.cosinefunc(6 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            float t = frameSize - 1;
            float[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
