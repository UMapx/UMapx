using System;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Grubb's statistic distribution.
    /// </summary>
    [Serializable]
    public class Grubb : IDistribution
    {
        #region Private data
        private GrubbDistribution distribution;
        #endregion

        #region Grubb components
        /// <summary>
        /// Initializes the Grubb's statistic distribution.
        /// </summary>
        /// <param name="samples">Number of samples.</param>
        public Grubb(int samples)
        {
            distribution = new GrubbDistribution(samples);
        }
        /// <summary>
        /// Gets the number of samples.
        /// </summary>
        public int Samples => distribution.NumberOfSamples;
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
        public float Mean => throw new NotSupportedException();
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance => throw new NotSupportedException();
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median => (float)distribution.Median;
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode => throw new NotSupportedException();
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
        public float Entropy => throw new NotSupportedException();
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x) => throw new NotSupportedException();
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x) => (float)distribution.DistributionFunction(x);
        #endregion
    }
}
