namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the interface for continuous wavelets.
    /// </summary>
    public interface IDoubleWavelet
    {
        #region Interface
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        double Scaling(double x);
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        double Wavelet(double x);
        #endregion
    }
}
