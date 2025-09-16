using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Tukey-Lambda distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Tukey_lambda_distribution
    /// </remarks>
    [Serializable]
    public class TukeyLambda : IDistribution
    {
        #region Private data
        private float lambda;
        #endregion

        #region Constructor
        /// <summary>
        /// Initializes the Tukey-Lambda distribution with the specified shape parameter.
        /// </summary>
        /// <param name="lambda">Shape parameter.</param>
        public TukeyLambda(float lambda = 0f)
        {
            this.lambda = lambda;
        }
        #endregion

        #region Distribution properties
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                if (lambda > 0f)
                {
                    return new RangeFloat(-1f / lambda, 1f / lambda);
                }

                return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean => 0f;
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (Maths.Abs(lambda) < 1e-6f)
                {
                    return (float)((Math.PI * Math.PI) / 3.0);
                }

                double l = lambda;
                double a = 2.0 / (l * l);
                double b = 1.0 / (1.0 + 2.0 * l);
                double c = Special.Gamma((float)(l + 1.0));
                double d = Special.Gamma((float)(2.0 * l + 2.0));
                double variance = a * (b - (c * c) / d);
                return (float)variance;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median => 0f;
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode => new float[] { 0f };
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness => 0f;
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        public float Excess => float.NaN;
        /// <summary>
        /// Gets the value of differential entropy.
        /// </summary>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion

        #region Methods
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Function(float x)
        {
            if (lambda > 0f)
            {
                float bound = 1f / lambda;
                if (x <= -bound || x >= bound)
                {
                    return 0f;
                }
            }

            double p = DistributionInternal(x);
            p = Maths.Range((float)p, 1e-10f, 1f - 1e-10f);
            double density = 1.0 / QuantileDensity(p);
            return (float)density;
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Distribution(float x)
        {
            if (lambda > 0f)
            {
                float bound = 1f / lambda;
                if (x <= -bound)
                {
                    return 0f;
                }

                if (x >= bound)
                {
                    return 1f;
                }
            }

            return (float)DistributionInternal(x);
        }
        #endregion

        #region Private helpers
        private double DistributionInternal(double x)
        {
            double lower = 0.0;
            double upper = 1.0;

            for (int i = 0; i < 64; i++)
            {
                double mid = 0.5 * (lower + upper);
                double value = Quantile(mid);

                if (value > x)
                {
                    upper = mid;
                }
                else
                {
                    lower = mid;
                }
            }

            return 0.5 * (lower + upper);
        }

        private double Quantile(double p)
        {
            if (p <= 0.0)
            {
                return lambda > 0f ? -1.0 / lambda : double.NegativeInfinity;
            }

            if (p >= 1.0)
            {
                return lambda > 0f ? 1.0 / lambda : double.PositiveInfinity;
            }

            if (Math.Abs(lambda) < 1e-12)
            {
                return Math.Log(p / (1.0 - p));
            }

            double a = Math.Pow(p, lambda);
            double b = Math.Pow(1.0 - p, lambda);
            return (a - b) / lambda;
        }

        private double QuantileDensity(double p)
        {
            if (Math.Abs(lambda) < 1e-12)
            {
                double numerator = 1.0;
                double denom = p * (1.0 - p);
                return numerator / denom;
            }

            double a = Math.Pow(p, lambda - 1.0);
            double b = Math.Pow(1.0 - p, lambda - 1.0);
            return a + b;
        }
        #endregion
    }
}
