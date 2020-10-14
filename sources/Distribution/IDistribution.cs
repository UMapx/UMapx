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
        RangeDouble Support { get; }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        double Mean { get; }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        double Variance { get; }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        double Median { get; }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        double Mode { get; }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        double Skewness { get; }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        double Excess { get; }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        double Entropy { get; }
        #endregion
    }
}
