using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Hann window function (Hanning).
    /// </summary>
    [Serializable]
    public class Hann : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Hann window function (Hanning).
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Hann(int frameSize)
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
            return (float)Math.Pow(Sine.sinefunc(x, frameSize), 2);
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
