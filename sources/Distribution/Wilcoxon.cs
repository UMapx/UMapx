using System;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Wilcoxon signed-rank distribution.
    /// </summary>
    [Serializable]
    public class Wilcoxon : IDistribution
    {
        #region Private data
        private WilcoxonDistribution distribution;
        #endregion

        #region Wilcoxon components
        /// <summary>
        /// Initializes the Wilcoxon distribution using the number of observations.
        /// </summary>
        /// <param name="samples">Number of observations.</param>
        public Wilcoxon(int samples)
        {
            distribution = new WilcoxonDistribution(samples);
        }
        /// <summary>
        /// Initializes the Wilcoxon distribution from rank statistics.
        /// </summary>
        /// <param name="ranks">Rank statistics.</param>
        /// <param name="exact">Whether to compute the exact distribution.</param>
        public Wilcoxon(float[] ranks, bool? exact = null)
        {
            distribution = new WilcoxonDistribution(DistributionHelper.ToDouble(ranks), exact);
        }
        /// <summary>
        /// Gets the number of effective samples.
        /// </summary>
        public int Samples => distribution.NumberOfSamples;
        /// <summary>
        /// Gets whether the distribution computes exact probabilities.
        /// </summary>
        public bool Exact => distribution.Exact;
        /// <summary>
        /// Gets or sets the continuity correction.
        /// </summary>
        public ContinuityCorrection Correction
        {
            get => distribution.Correction;
            set => distribution.Correction = value;
        }
        /// <summary>
        /// Gets the statistic table for the exact distribution.
        /// </summary>
        public float[] Table => DistributionHelper.ToFloat(distribution.Table);
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
