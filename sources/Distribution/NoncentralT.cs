using System;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the noncentral t-distribution.
    /// </summary>
    [Serializable]
    public class NoncentralT : IDistribution
    {
        #region Private data
        private NoncentralTDistribution distribution;
        #endregion

        #region NoncentralT components
        /// <summary>
        /// Initializes the noncentral t-distribution.
        /// </summary>
        /// <param name="degreesOfFreedom">Degrees of freedom.</param>
        /// <param name="noncentrality">Noncentrality parameter.</param>
        public NoncentralT(float degreesOfFreedom, float noncentrality)
        {
            distribution = new NoncentralTDistribution(degreesOfFreedom, noncentrality);
        }
        /// <summary>
        /// Gets the degrees of freedom.
        /// </summary>
        public float DegreesOfFreedom => (float)distribution.DegreesOfFreedom;
        /// <summary>
        /// Gets the noncentrality parameter.
        /// </summary>
        public float Noncentrality => (float)distribution.Noncentrality;
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
        /// Gets the mean value.
        /// </summary>
        public float Mean => (float)distribution.Mean;
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance => (float)distribution.Variance;
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median => (float)distribution.Median;
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode => new[] { (float)distribution.Mode };
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
