using System;
using UMapx.Core;

namespace UMapx.Window
{
    /// <summary>
    /// Defines a generalized window normal function.
    /// </summary>
    [Serializable]
    public class Normal : WindowBase
    {
        #region Private data
        private float sigma = 1;
        private float p = 2;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes a generalized window normal function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="sigma">Standard deviation (>0)</param>
        /// <param name="pow">Power<remarks>For p = 2 - Gaussian window</remarks></param>
        public Normal(int frameSize, float sigma = 1, float pow = 2)
        {
            this.Sigma = sigma;
            this.FrameSize = frameSize;
            this.p = pow;
        }
        /// <summary>
        /// Gets or sets the standard deviation (>0).
        /// </summary>
        public float Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Power.
        /// </summary>
        public float Pow
        {
            get
            {
                return this.p;
            }
            set
            {
                this.p = value;
            }
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Value</returns>
        public override float Function(float x, int frameSize)
        {
            float a = (frameSize - 1) / 2;
            float t = (x - a) / (sigma * a);
            return Maths.Exp(-Maths.Pow(Maths.Abs(t), p));
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
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
