using System;
using Accord.Statistics.Distributions.Univariate;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the skew normal distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Skew_normal_distribution
    /// </remarks>
    [Serializable]
    public class SkewNormal : IDistribution
    {
        #region Private data
        private SkewNormalDistribution distribution;
        #endregion

        #region SkewNormal components
        /// <summary>
        /// Initializes the skew normal distribution.
        /// </summary>
        public SkewNormal()
        {
            distribution = new SkewNormalDistribution();
        }
        /// <summary>
        /// Initializes the skew normal distribution.
        /// </summary>
        /// <param name="location">Location parameter.</param>
        public SkewNormal(float location)
        {
            distribution = new SkewNormalDistribution(location);
        }
        /// <summary>
        /// Initializes the skew normal distribution.
        /// </summary>
        /// <param name="location">Location parameter.</param>
        /// <param name="scale">Scale parameter (greater than zero).</param>
        public SkewNormal(float location, float scale)
        {
            distribution = new SkewNormalDistribution(location, scale);
        }
        /// <summary>
        /// Initializes the skew normal distribution.
        /// </summary>
        /// <param name="location">Location parameter.</param>
        /// <param name="scale">Scale parameter (greater than zero).</param>
        /// <param name="shape">Shape parameter.</param>
        public SkewNormal(float location, float scale, float shape)
        {
            distribution = new SkewNormalDistribution(location, scale, shape);
        }
        /// <summary>
        /// Gets or sets the location parameter.
        /// </summary>
        public float Location
        {
            get => (float)distribution.Location;
            set => distribution = new SkewNormalDistribution(value, Scale, Shape);
        }
        /// <summary>
        /// Gets or sets the scale parameter.
        /// </summary>
        public float Scale
        {
            get => (float)distribution.Scale;
            set => distribution = new SkewNormalDistribution(Location, value, Shape);
        }
        /// <summary>
        /// Gets or sets the shape parameter.
        /// </summary>
        public float Shape
        {
            get => (float)distribution.Shape;
            set => distribution = new SkewNormalDistribution(Location, Scale, value);
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
        public float Skewness => (float)distribution.Skewness;
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
                double alpha = Shape;
                double delta = alpha / Math.Sqrt(1.0 + alpha * alpha);
                double b = delta * Math.Sqrt(2.0 / Math.PI);
                double c = 1.0 - 2.0 * delta * delta / Math.PI;
                double excess = 2.0 * (Math.PI - 3.0) * b * b * b * b / (c * c);
                return (float)excess;
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
