using System;
using Accord.Math;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines a general continuous distribution based on provided density or distribution functions.
    /// </summary>
    [Serializable]
    public class GeneralContinuous : IDistribution
    {
        #region Private data
        private GeneralContinuousDistribution distribution;
        #endregion

        #region GeneralContinuous components
        /// <summary>
        /// Initializes the distribution from a density function defined over the entire real line.
        /// </summary>
        /// <param name="density">Density function.</param>
        public GeneralContinuous(Func<double, double> density)
        {
            distribution = GeneralContinuousDistribution.FromDensityFunction(density);
        }
        /// <summary>
        /// Initializes the distribution from a density function over a support range.
        /// </summary>
        /// <param name="support">Support range.</param>
        /// <param name="density">Density function.</param>
        public GeneralContinuous(RangeFloat support, Func<double, double> density)
        {
            distribution = GeneralContinuousDistribution.FromDensityFunction(new DoubleRange(support.Min, support.Max), density);
        }
        /// <summary>
        /// Initializes the distribution from density and distribution functions.
        /// </summary>
        /// <param name="support">Support range.</param>
        /// <param name="density">Density function.</param>
        /// <param name="distributionFunction">Cumulative distribution function.</param>
        public GeneralContinuous(RangeFloat support, Func<double, double> density, Func<double, double> distributionFunction)
        {
            distribution = new GeneralContinuousDistribution(new DoubleRange(support.Min, support.Max), density, distributionFunction);
        }
        /// <summary>
        /// Initializes the distribution from an existing Accord distribution.
        /// </summary>
        /// <param name="distribution">Accord distribution.</param>
        public GeneralContinuous(GeneralContinuousDistribution distribution)
        {
            this.distribution = distribution ?? throw new ArgumentNullException(nameof(distribution));
        }
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
