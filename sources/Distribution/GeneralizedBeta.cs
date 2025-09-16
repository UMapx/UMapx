using System;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the generalized beta distribution.
    /// </summary>
    [Serializable]
    public class GeneralizedBeta : IDistribution
    {
        #region Private data
        private GeneralizedBetaDistribution distribution;
        #endregion

        #region GeneralizedBeta components
        /// <summary>
        /// Initializes the generalized beta distribution defined on (0, 1).
        /// </summary>
        /// <param name="alpha">Shape parameter α.</param>
        /// <param name="beta">Shape parameter β.</param>
        public GeneralizedBeta(float alpha, float beta)
        {
            distribution = new GeneralizedBetaDistribution(alpha, beta);
        }
        /// <summary>
        /// Initializes the generalized beta distribution on (min, max).
        /// </summary>
        /// <param name="alpha">Shape parameter α.</param>
        /// <param name="beta">Shape parameter β.</param>
        /// <param name="min">Minimum value.</param>
        /// <param name="max">Maximum value.</param>
        public GeneralizedBeta(float alpha, float beta, float min, float max)
        {
            distribution = new GeneralizedBetaDistribution(alpha, beta, min, max);
        }
        /// <summary>
        /// Initializes the beta PERT distribution.
        /// </summary>
        /// <param name="min">Minimum value.</param>
        /// <param name="max">Maximum value.</param>
        /// <param name="mode">Most probable value.</param>
        public GeneralizedBeta(float min, float max, float mode)
        {
            distribution = GeneralizedBetaDistribution.Pert(min, max, mode);
        }
        /// <summary>
        /// Initializes the beta PERT distribution with scale.
        /// </summary>
        /// <param name="min">Minimum value.</param>
        /// <param name="max">Maximum value.</param>
        /// <param name="mode">Most probable value.</param>
        /// <param name="scale">Scale (λ) parameter.</param>
        public GeneralizedBeta(float min, float max, float mode, float scale)
        {
            distribution = GeneralizedBetaDistribution.Pert(min, max, mode, scale);
        }
        /// <summary>
        /// Gets the shape parameter α.
        /// </summary>
        public float Alpha => (float)distribution.Alpha;
        /// <summary>
        /// Gets the shape parameter β.
        /// </summary>
        public float Beta => (float)distribution.Beta;
        /// <summary>
        /// Gets the minimum value.
        /// </summary>
        public float Min => (float)distribution.Min;
        /// <summary>
        /// Gets the maximum value.
        /// </summary>
        public float Max => (float)distribution.Max;
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support => new RangeFloat(Min, Max);
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
        public float Skewness
        {
            get
            {
                double a = Alpha;
                double b = Beta;
                double numerator = 2.0 * (b - a) * Math.Sqrt(a + b + 1.0);
                double denominator = (a + b + 2.0) * Math.Sqrt(a * b);
                return (float)(numerator / denominator);
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get
            {
                double a = Alpha;
                double b = Beta;
                double numerator = 6.0 * ((a - b) * (a - b) * (a + b + 1.0) - a * b * (a + b + 2.0));
                double denominator = a * b * (a + b + 2.0) * (a + b + 3.0);
                return (float)(numerator / denominator);
            }
        }
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
