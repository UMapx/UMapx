using System;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous Mexican hat wavelet.
    /// </summary>
    [Serializable]
    public class MexicanHatWavelet : IDoubleWavelet
    {
        #region Wavelet components
        /// <summary>
        /// Initializes the continuous Mexican hat wavelet.
        /// </summary>
        public MexicanHatWavelet() { }
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
            return 2.0 / (Math.Sqrt(3) * Math.Pow(Math.PI, 0.25)) * (1 - x2) * Math.Exp(-x2 / 2);
        }
        #endregion
    }
}
