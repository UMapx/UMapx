using System;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the empirical distribution.
    /// </summary>
    [Serializable]
    public class Empirical : IDistribution
    {
        #region Private data
        private EmpiricalDistribution distribution;
        #endregion

        #region Empirical components
        /// <summary>
        /// Initializes the empirical distribution from samples.
        /// </summary>
        /// <param name="samples">Samples.</param>
        public Empirical(float[] samples)
        {
            distribution = new EmpiricalDistribution(DistributionHelper.ToDouble(samples));
        }
        /// <summary>
        /// Initializes the empirical distribution from samples with smoothing parameter.
        /// </summary>
        /// <param name="samples">Samples.</param>
        /// <param name="smoothing">Smoothing parameter.</param>
        public Empirical(float[] samples, float smoothing)
        {
            distribution = new EmpiricalDistribution(DistributionHelper.ToDouble(samples), smoothing);
        }
        /// <summary>
        /// Initializes the empirical distribution from samples with fractional weights.
        /// </summary>
        /// <param name="samples">Samples.</param>
        /// <param name="weights">Fractional weights.</param>
        public Empirical(float[] samples, float[] weights)
        {
            distribution = new EmpiricalDistribution(DistributionHelper.ToDouble(samples), DistributionHelper.ToDouble(weights));
        }
        /// <summary>
        /// Initializes the empirical distribution from samples with repetition weights.
        /// </summary>
        /// <param name="samples">Samples.</param>
        /// <param name="weights">Repetition weights.</param>
        public Empirical(float[] samples, int[] weights)
        {
            distribution = new EmpiricalDistribution(DistributionHelper.ToDouble(samples), weights);
        }
        /// <summary>
        /// Initializes the empirical distribution from samples with fractional weights and smoothing parameter.
        /// </summary>
        /// <param name="samples">Samples.</param>
        /// <param name="weights">Fractional weights.</param>
        /// <param name="smoothing">Smoothing parameter.</param>
        public Empirical(float[] samples, float[] weights, float smoothing)
        {
            distribution = new EmpiricalDistribution(DistributionHelper.ToDouble(samples), DistributionHelper.ToDouble(weights), smoothing);
        }
        /// <summary>
        /// Initializes the empirical distribution from samples with repetition weights and smoothing parameter.
        /// </summary>
        /// <param name="samples">Samples.</param>
        /// <param name="weights">Repetition weights.</param>
        /// <param name="smoothing">Smoothing parameter.</param>
        public Empirical(float[] samples, int[] weights, float smoothing)
        {
            distribution = new EmpiricalDistribution(DistributionHelper.ToDouble(samples), weights, smoothing);
        }
        /// <summary>
        /// Gets smoothing parameter.
        /// </summary>
        public float Smoothing => (float)distribution.Smoothing;
        /// <summary>
        /// Gets samples.
        /// </summary>
        public float[] Samples => DistributionHelper.ToFloat(distribution.Samples);
        /// <summary>
        /// Gets fractional weights.
        /// </summary>
        public float[] Weights => DistributionHelper.ToFloat(distribution.Weights);
        /// <summary>
        /// Gets repetition counts.
        /// </summary>
        public int[] Counts => distribution.Counts;
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
