using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Shannon wavelet.
    /// </summary>
    [Serializable]
    public class ShannonWavelet : IDoubleWavelet
    {
        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Shannon wavelet.
        /// </summary>
        public ShannonWavelet() { }
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
            double t = x / 2;
            return Special.Sinc(t) * Math.Cos(3 * Math.PI * t);
        }
        #endregion
    }
}
