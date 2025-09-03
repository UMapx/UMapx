using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Parzen window function.
    /// </summary>
    [Serializable]
    public class Parzen : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Parzen window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Parzen(int frameSize)
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
            // coefficients:
            float y = Math.Abs(x);
            float c = frameSize / 2.0f;
            float a = x / c;
            float b = y / c;

            // props:
            if ((y >= 0) &&
                (y <= frameSize / 4))
            {
                return 1.0f - 6.0f * a * a * (1.0f - b);
            }
            else if (y >= frameSize / 4 &&
                y <= frameSize / 2)
            {
                return 2 * (float)Math.Pow(1 - b, 3);
            }
            return 0.0f;
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override float[] GetWindow(int frameSize)
        {
            float t = (frameSize - 1) / 2.0f;
            float[] x = MatrixF.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
}
