using System;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Mann-Whitney U statistic distribution.
    /// </summary>
    [Serializable]
    public class MannWhitney : IDistribution
    {
        #region Private data
        private MannWhitneyDistribution distribution;
        #endregion

        #region MannWhitney components
        /// <summary>
        /// Initializes the distribution using sample sizes.
        /// </summary>
        /// <param name="n1">Number of observations in the first sample.</param>
        /// <param name="n2">Number of observations in the second sample.</param>
        public MannWhitney(int n1, int n2)
        {
            distribution = new MannWhitneyDistribution(n1, n2);
        }
        /// <summary>
        /// Initializes the distribution from rank statistics.
        /// </summary>
        /// <param name="ranks">Global rank statistics.</param>
        /// <param name="n1">Number of observations in the first sample.</param>
        /// <param name="n2">Number of observations in the second sample.</param>
        /// <param name="exact">Whether to compute the exact distribution.</param>
        public MannWhitney(float[] ranks, int n1, int n2, bool? exact = null)
        {
            distribution = new MannWhitneyDistribution(DistributionHelper.ToDouble(ranks), n1, n2, exact);
        }
        /// <summary>
        /// Initializes the distribution from rank statistics for both samples.
        /// </summary>
        /// <param name="ranks1">Ranks for the first sample.</param>
        /// <param name="ranks2">Ranks for the second sample.</param>
        /// <param name="exact">Whether to compute the exact distribution.</param>
        public MannWhitney(float[] ranks1, float[] ranks2, bool? exact = null)
        {
            distribution = new MannWhitneyDistribution(DistributionHelper.ToDouble(ranks1), DistributionHelper.ToDouble(ranks2), exact);
        }
        /// <summary>
        /// Gets the number of observations in the first sample.
        /// </summary>
        public int Samples1 => distribution.NumberOfSamples1;
        /// <summary>
        /// Gets the number of observations in the second sample.
        /// </summary>
        public int Samples2 => distribution.NumberOfSamples2;
        /// <summary>
        /// Gets or sets the continuity correction mode.
        /// </summary>
        public ContinuityCorrection Correction
        {
            get => distribution.Correction;
            set => distribution.Correction = value;
        }
        /// <summary>
        /// Gets whether this distribution computes exact probabilities.
        /// </summary>
        public bool Exact => distribution.Exact;
        /// <summary>
        /// Gets the statistic table used for the exact distribution.
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
