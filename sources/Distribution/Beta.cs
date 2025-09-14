using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the beta distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Beta_distribution
    /// </remarks>
    [Serializable]
    public class Beta : IDistribution
    {
        #region Private data
        private float a = 1;
        private float b = 1;
        #endregion

        #region Beta components
        /// <summary>
        /// Initializes beta distribution.
        /// </summary>
        public Beta() { }
        /// <summary>
        /// Initializes beta distribution.
        /// </summary>
        /// <param name="a">Parameter a</param>
        /// <param name="b">Parameter b</param>
        public Beta(float a, float b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the parameter a.
        /// </summary>
        public float A
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Parameter a must be greater than zero");
                }

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the parameter b.
        /// </summary>
        public float B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentException("Parameter b must be greater than zero");
                }

                this.b = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(0, 1);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return a / (a + b);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return (a * b) / Maths.Pow(a + b) / (a + b + 1);
            }
        }
        /// <summary>
        /// Gets the mode value. If <c>a ≤ 1</c> and <c>b > 1</c>,
        /// the mode is at 0. If <c>b ≤ 1</c> and <c>a > 1</c>,
        /// the mode is at 1. When both parameters are greater than
        /// one, the mode is calculated as <c>(a - 1) / (a + b - 2)</c>;
        /// otherwise, the mode is undefined and returns <see cref="float.NaN" />.
        /// </summary>
        public float Mode
        {
            get
            {
                if (a <= 1 && b > 1)
                {
                    return 0f;
                }
                if (b <= 1 && a > 1)
                {
                    return 1f;
                }
                if (a > 1 && b > 1)
                {
                    return (a - 1) / (a + b - 2);
                }
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return 2 * (b - a) * Maths.Sqrt(a + b + 1) / (a + b + 2) / Maths.Sqrt(a * b);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                float a2 = a * a, b2 = b * b, a3 = a2 * a;
                return 6 * (a3 - a2 * (2 * b - 1) + b2 * (b + 1) - 2 * a * b * (b + 2)) / (a * b * (a + b + 2) * (a + b + 3));
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x > 1)
            {
                return 0;
            }
            else if (x < 0)
            {
                return 0;
            }
            return Maths.Pow(x, a - 1) * Maths.Pow(1 - x, b - 1) / Special.Beta(a, b);
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// Uses the regularized incomplete beta function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        /// <example>
        /// Beta beta = new Beta(2f, 3f);
        /// float cdf = beta.Distribution(0.5f);
        /// </example>
        public float Distribution(float x)
        {
            if (x > 1)
            {
                return 1;
            }
            else if (x < 0)
            {
                return 0;
            }
            return Special.BetaIncompleteRegularized(a, b, x);
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
    }
}
