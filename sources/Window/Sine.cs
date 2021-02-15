using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the sine window function.
    /// </summary>
    [Serializable]
    public class Sine : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the sine window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Sine(int frameSize)
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
            return (float)Math.Sin(Math.PI * x / (frameSize - 1));
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

        #region Static components
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>float precision floating point number</returns>
        internal static float sinefunc(float x, int frameSize)
        {
            return (float)Math.Sin(Math.PI * x / (frameSize - 1));
        }
        #endregion
    }
}
