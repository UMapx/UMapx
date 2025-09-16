using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the folded normal distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Folded_normal_distribution
    /// </remarks>
    [Serializable]
    public class FoldedNormal : IDistribution
    {
        #region Private data
        private float mu = 0f;
        private float sigma = 1f;
        #endregion

        #region Folded normal components
        /// <summary>
        /// Initializes the folded normal distribution.
        /// </summary>
        public FoldedNormal() { }
        /// <summary>
        /// Initializes the folded normal distribution.
        /// </summary>
        /// <param name="mu">Mean value of the underlying Gaussian variable.</param>
        /// <param name="sigma">Standard deviation of the underlying Gaussian variable.</param>
        public FoldedNormal(float mu, float sigma)
        {
            Mu = mu;
            Sigma = sigma;
        }
        /// <summary>
        /// Gets or sets the mean value of the underlying Gaussian variable.
        /// </summary>
        public float Mu
        {
            get { return this.mu; }
            set { this.mu = value; }
        }
        /// <summary>
        /// Gets or sets the standard deviation of the underlying Gaussian variable.
        /// </summary>
        public float Sigma
        {
            get { return this.sigma; }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(0, float.PositiveInfinity); }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                double sigma2 = sigma * sigma;
                double coeff = sigma * Math.Sqrt(2.0 / Math.PI) * Math.Exp(-(mu * mu) / (2.0 * sigma2));
                double corr = mu * (1.0 - 2.0 * StandardNormalCdf(-mu / sigma));
                return (float)(coeff + corr);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                double mean = Mean;
                return (float)(mu * mu + sigma * sigma - mean * mean);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return FindQuantile(0.5f);
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get { return new float[] { ComputeMode() }; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { return float.NaN; }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get { return float.NaN; }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < 0)
            {
                return 0f;
            }

            double a = (x - mu) / sigma;
            double b = (-x - mu) / sigma;
            return (float)((StandardNormalPdf(a) + StandardNormalPdf(b)) / sigma);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= 0)
            {
                return 0f;
            }

            double den = Math.Sqrt(2.0) * sigma;
            double a = (x + mu) / den;
            double b = (x - mu) / den;
            return (float)(0.5 * (Special.Erf((float)a) + Special.Erf((float)b)));
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
        private static double StandardNormalPdf(double z)
        {
            return Math.Exp(-0.5 * z * z) / Math.Sqrt(2.0 * Math.PI);
        }

        private static double StandardNormalCdf(double z)
        {
            return 0.5 * (1.0 + Special.Erf((float)(z / Math.Sqrt(2.0))));
        }

        private float ComputeMode()
        {
            if (sigma == 0)
            {
                return 0f;
            }

            double absMu = Math.Abs(mu);
            if (absMu == 0)
            {
                return 0f;
            }

            double ratio = absMu / sigma;
            double exponent = -absMu * absMu / (2.0 * sigma * sigma);
            double arg = ratio * ratio * Math.Exp(exponent);
            if (arg <= 0)
            {
                return 0f;
            }

            double w = Special.LambertW((float)arg);
            double value = sigma / Math.Sqrt(2.0) * Math.Sqrt(w + ratio * ratio);
            return (float)value;
        }

        private float FindQuantile(float p)
        {
            if (p <= 0f)
                return 0f;
            if (p >= 1f)
                return float.PositiveInfinity;

            double lower = 0.0;
            double upper = Math.Abs(mu) + 10.0 * sigma + 10.0;

            for (int i = 0; i < 50; i++)
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
