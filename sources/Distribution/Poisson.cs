using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Poisson distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// <see href="https://en.wikipedia.org/wiki/Poisson_distribution"/>.
    /// </remarks>
    [Serializable]
    public class Poisson : IDistribution
    {
        #region Private data
        private float l = 1;
        #endregion

        #region Poisson components
        /// <summary>
        /// Initializes the Poisson distribution.
        /// </summary>
        public Poisson() { }
        /// <summary>
        /// Initializes the Poisson distribution.
        /// </summary>
        /// <param name="lambda">Parameter λ (0, +inf).</param>
        public Poisson(float lambda)
        {
            Lambda = lambda;
        }
        /// <summary>
        /// Gets or sets the value of the parameter λ (0, +inf).
        /// </summary>
        public float Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Lambda must be positive");

                this.l = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(0, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return l;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return l;
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                float mode = Maths.Floor(l);

                if (mode > 0f && l == mode)
                {
                    return new float[] { mode - 1f, mode };
                }

                return new float[] { mode };
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        /// <remarks>
        /// Returns the lower median by cumulative-probability search, rounded to binary32.
        /// </remarks>
        public float Median
        {
            get
            {
                // Beyond this limit binary32 cannot resolve the sub-unit median correction.
                if (l >= 16777216) return l;
                int low = 0, high = (int)Math.Ceiling((double)l + 1);
                while (low < high)
                {
                    int mid = low + (high - low) / 2;
                    if (Special.DistributionGamma(mid + 1.0, l, true) >= 0.5) high = mid;
                    else low = mid + 1;
                }
                return low;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return Maths.Pow(l, -0.5f);
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
                return Maths.Pow(l, -1.0f);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Function(float x)
        {
            if (float.IsNaN(x)) return float.NaN;
            if (x < 0 || float.IsPositiveInfinity(x) || x != Math.Floor(x)) return 0;
            return (float)Math.Exp(DistributionNumerics.PoissonLogMass(l, x));
        }
        /// <summary>
        /// Returns the value of the probability mass cumulative function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Distribution(float x)
        {
            if (float.IsNaN(x)) return float.NaN;
            if (x < 0) return 0;
            if (float.IsPositiveInfinity(x)) return 1;
            return (float)Special.DistributionGamma(Math.Floor(x) + 1, l, true);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value.</returns>
        public float Entropy
        {
            get
            {
                double lambda = l;
                if (lambda > 1000)
                {
                    double inverse = 1 / lambda;
                    return (float)(0.5 + DistributionNumerics.LogSqrtTwoPi + 0.5 * Math.Log(lambda)
                        - inverse / 12 - inverse * inverse / 24 - 19 * inverse * inverse * inverse / 360);
                }
                int mode = (int)Math.Floor(lambda);
                double logMass = DistributionNumerics.PoissonLogMass(lambda, mode);
                double atMode = Math.Exp(logMass), sum = -atMode * logMass;
                double mass = atMode;
                for (int k = mode; k > 0; k--)
                {
                    mass *= k / lambda;
                    if (mass > 0) sum -= mass * Math.Log(mass);
                }
                mass = atMode;
                int limit = mode + (int)Math.Ceiling(14 * Math.Sqrt(lambda)) + 50;
                for (int k = mode + 1; k <= limit; k++)
                {
                    mass *= lambda / k;
                    if (mass > 0) sum -= mass * Math.Log(mass);
                }
                return (float)sum;
            }
        }
        #endregion
    }
}
