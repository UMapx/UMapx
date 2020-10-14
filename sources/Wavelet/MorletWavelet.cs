using System;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Morlet wavelet.
    /// </summary>
    [Serializable]
    public class MorletWavelet : IDoubleWavelet
    {
        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Morlet wavelet.
        /// </summary>
        public MorletWavelet()
        {
        }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public double Scaling(double x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public double Wavelet(double x)
        {
            double x2 = x * x;
            return Math.Exp(-x2 / 2) * Math.Cos(5 * x);
        }
        #endregion
    }
}
