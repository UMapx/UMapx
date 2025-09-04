using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the cosine window function.
    /// </summary>
    [Serializable]
    public class Cosine : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the cosine window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Cosine(int frameSize)
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
            return (float)Math.Cos(Maths.Pi * x / (frameSize - 1));
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            float t = (frameSize - 1) / 2.0f;
            float[] x = Matrix.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion

        #region Static components
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Factor</returns>
        internal static float Cosinefunc(float x, int frameSize)
        {
            return (float)Math.Cos(Maths.Pi * x / (frameSize - 1));
        }
        #endregion
    }
}
