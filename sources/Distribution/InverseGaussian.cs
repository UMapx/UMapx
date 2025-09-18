using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the inverse Gaussian (Wald) distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
    /// </remarks>
    [Serializable]
    public class InverseGaussian : IDistribution
    {
        #region Private data
        private float mu = 1f;
        private float lambda = 1f;
        #endregion

        #region Inverse Gaussian components
        /// <summary>
        /// Initializes the inverse Gaussian distribution.
        /// </summary>
        public InverseGaussian() { }
        /// <summary>
        /// Initializes the inverse Gaussian distribution.
        /// </summary>
        /// <param name="mu">Mean parameter μ (0, +inf)</param>
        /// <param name="lambda">Shape parameter λ (0, +inf)</param>
        public InverseGaussian(float mu, float lambda)
        {
            Mu = mu;
            Lambda = lambda;
        }
        /// <summary>
        /// Gets or sets the mean parameter μ (0, +inf).
        /// </summary>
        public float Mu
        {
            get { return this.mu; }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the shape parameter λ (0, +inf).
        /// </summary>
        public float Lambda
        {
            get { return this.lambda; }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.lambda = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(0f, float.PositiveInfinity); }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return this.mu; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return mu * mu * mu / lambda; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { return FindQuantile(0.5f); }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                double ratio = (3.0 * mu) / (2.0 * lambda);
                double term = Math.Sqrt(1.0 + (9.0 * mu * mu) / (4.0 * lambda * lambda));
                double mode = mu * (term - ratio);
                return new float[] { (float)mode };
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { return 3f * Maths.Sqrt(mu / lambda); }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get { return 15f * mu / lambda; }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= 0f)
                return 0f;

            double coef = Math.Sqrt(lambda / (2.0 * Math.PI * x * x * x));
            double expo = -lambda * Math.Pow(x - mu, 2) / (2.0 * mu * mu * x);
            return (float)(coef * Math.Exp(expo));
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= 0f)
                return 0f;

            double sqrt = Math.Sqrt(lambda / x);
            double a = 0.5 * Special.Erfc((float)(sqrt * (mu - x) / (Maths.Sqrt2 * mu)));
            double b = 0.5 * Special.Erfc((float)(sqrt * (mu + x) / (Maths.Sqrt2 * mu)));
            double c = Math.Exp((2.0 * lambda) / mu);
            return (float)(a + b * c);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion

        #region Private methods
        private float FindQuantile(float p)
        {
            if (p <= 0f)
                return 0f;
            if (p >= 1f)
                return float.PositiveInfinity;

            double lower = 0.0;
            double upper = mu + 10.0 * Math.Sqrt(Variance) + 10.0;

            for (int i = 0; i < 60; i++)
            {
                double mid = 0.5 * (lower + upper);
                double cdf = Distribution((float)mid);
                if (cdf > p)
                    upper = mid;
                else
                    lower = mid;
            }

            return (float)(0.5 * (lower + upper));
        }
        #endregion
    }
}
