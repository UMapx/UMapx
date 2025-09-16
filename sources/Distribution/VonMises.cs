using System;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the von-Mises (circular normal) distribution.
    /// </summary>
    [Serializable]
    public class VonMises : IDistribution
    {
        #region Private data
        private VonMisesDistribution distribution;
        #endregion

        #region VonMises components
        /// <summary>
        /// Initializes the von-Mises distribution with zero mean and unit concentration.
        /// </summary>
        public VonMises()
        {
            distribution = new VonMisesDistribution();
        }
        /// <summary>
        /// Initializes the von-Mises distribution with zero mean.
        /// </summary>
        /// <param name="concentration">Concentration parameter κ.</param>
        public VonMises(float concentration)
        {
            distribution = new VonMisesDistribution(concentration);
        }
        /// <summary>
        /// Initializes the von-Mises distribution.
        /// </summary>
        /// <param name="mean">Mean angle μ.</param>
        /// <param name="concentration">Concentration parameter κ.</param>
        public VonMises(float mean, float concentration)
        {
            distribution = new VonMisesDistribution(mean, concentration);
        }
        /// <summary>
        /// Gets the mean value μ.
        /// </summary>
        public float Mean => (float)distribution.Mean;
        /// <summary>
        /// Gets the median value μ.
        /// </summary>
        public float Median => (float)distribution.Median;
        /// <summary>
        /// Gets the mode value μ.
        /// </summary>
        public float[] Mode => new[] { (float)distribution.Mode };
        /// <summary>
        /// Gets the concentration parameter κ.
        /// </summary>
        public float Concentration => (float)distribution.Concentration;
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                var support = distribution.Support;
                return new RangeFloat((float)support.Min, (float)support.Max);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance => (float)distribution.Variance;
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness => throw new NotSupportedException();
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess => throw new NotSupportedException();
        /// <summary>
        /// Gets the value of differential entropy.
        /// </summary>
        public float Entropy => (float)distribution.Entropy;
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x) => (float)distribution.ProbabilityDensityFunction(x);
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x) => (float)distribution.DistributionFunction(x);
        #endregion
    }
}
