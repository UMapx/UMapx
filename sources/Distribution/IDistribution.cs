using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the distribution interface.
    /// </summary>
    public interface IDistribution
    {
        #region Components
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        RangeFloat Support { get; }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        float Mean { get; }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        float Variance { get; }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        float Median { get; }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        float Mode { get; }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        float Skewness { get; }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        float Excess { get; }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        float Entropy { get; }
        #endregion
    }
}
