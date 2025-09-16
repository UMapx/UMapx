using System;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the empirical hazard distribution.
    /// </summary>
    [Serializable]
    public class EmpiricalHazard : IDistribution
    {
        #region Private data
        private EmpiricalHazardDistribution distribution;
        #endregion

        #region EmpiricalHazard components
        /// <summary>
        /// Initializes the empirical hazard distribution.
        /// </summary>
        public EmpiricalHazard()
        {
            distribution = new EmpiricalHazardDistribution();
        }
        /// <summary>
        /// Initializes the empirical hazard distribution with times and hazard values.
        /// </summary>
        /// <param name="times">Times.</param>
        /// <param name="hazards">Hazard rates.</param>
        public EmpiricalHazard(float[] times, float[] hazards)
        {
            distribution = new EmpiricalHazardDistribution(DistributionHelper.ToDouble(times), DistributionHelper.ToDouble(hazards));
        }
        /// <summary>
        /// Initializes the empirical hazard distribution with estimator.
        /// </summary>
        /// <param name="times">Times.</param>
        /// <param name="hazards">Hazard rates.</param>
        /// <param name="estimator">Survival estimator.</param>
        public EmpiricalHazard(float[] times, float[] hazards, SurvivalEstimator estimator)
        {
            distribution = new EmpiricalHazardDistribution(DistributionHelper.ToDouble(times), DistributionHelper.ToDouble(hazards), estimator);
        }
        /// <summary>
        /// Initializes the empirical hazard distribution with estimator.
        /// </summary>
        /// <param name="estimator">Survival estimator.</param>
        public EmpiricalHazard(SurvivalEstimator estimator)
        {
            distribution = new EmpiricalHazardDistribution(estimator);
        }
        /// <summary>
        /// Gets the survival estimator.
        /// </summary>
        public SurvivalEstimator Estimator => distribution.Estimator;
        /// <summary>
        /// Gets the times.
        /// </summary>
        public float[] Times => DistributionHelper.ToFloat(distribution.Times);
        /// <summary>
        /// Gets the hazard values.
        /// </summary>
        public float[] Hazards => DistributionHelper.ToFloat(distribution.Hazards);
        /// <summary>
        /// Gets the survival function values.
        /// </summary>
        public float[] Survivals => DistributionHelper.ToFloat(distribution.Survivals);
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
