using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Poisson wavelet.
    /// </summary>
    [Serializable]
    public class PoissonWavelet : IFloatWavelet
    {
        #region Private data
        private int n;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Poisson wavelet.
        /// </summary>
        /// <param name="n">Order [1, +inf)</param>
        public PoissonWavelet(int n = 1)
        {
            N = n;
        }
        /// <summary>
        /// Gets or sets the order [1, +inf).
        /// </summary>
        public int N
        {
            get
            {
                return n;
            }
            set
            {
                if (value < 1)
                    throw new Exception("Invalid argument value");

                n = value;
            }
        }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public float Scaling(float x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public float Wavelet(float x)
        {
            if (x < 0)
            {
                return 0;
            }
            return (x - n) / Special.Factorial(n) * (float)Math.Pow(x, n - 1) * (float)Math.Exp(-x);
        }
        #endregion
    }
}
