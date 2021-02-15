namespace UMapx.Response
{
    /// <summary>
    /// Defines the general interface of response Filters.
    /// </summary>
    public interface IResponse
    {
        #region Interface Components
        /// <summary>
        /// Returns an array of filter response values when a discrete function is supplied.
        /// </summary>
        /// <param name="u">Array</param>
        /// <returns>Discrete function in a Cartesian coordinate system</returns>
        float[] Reaction(float[] u);
        /// <summary>
        /// Returns the frequency response of the filter.
        /// </summary>
        /// <param name="w">Array of frequencies (rad / s)</param>
        /// <returns>Discrete function in a Cartesian coordinate system</returns>
        float[] Amplitude(float[] w);
        /// <summary>
        /// Returns the phase-frequency response of a filter.
        /// </summary>
        /// <param name="w">Array of frequencies (rad / s)</param>
        /// <returns>Discrete function in a Cartesian coordinate system</returns>
        float[] Phase(float[] w);
        /// <summary>
        /// Returns the amplitude value at the given frequency.
        /// </summary>
        /// <param name="w">Frequency (rad / s)</param>
        /// <returns>float precision floating point number</returns>
        float Amplitude(float w);
        /// <summary>
        /// Returns the phase value at the given frequency.
        /// </summary>
        /// <param name="w">Frequency (rad / s)</param>
        /// <returns>float precision floating point number</returns>
        float Phase(float w);
        #endregion
    }
}
