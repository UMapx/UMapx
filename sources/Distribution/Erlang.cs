using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the distribution of Erlang.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Erlang_distribution
    /// </remarks>
    [Serializable]
    public class Erlang : IDistribution
    {
        #region Private data
        private int k = 1;
        private float lambda = 0.5f;
        #endregion

        #region Erlang distribution
        /// <summary>
        /// Initializes the distribution of Erlang.
        /// </summary>
        /// <param name="k">Form parameter k ∈ (0, +inf)</param>
        /// <param name="lambda">λ-parameter λ ∈ (0, +inf)</param>
        public Erlang(int k, float lambda)
        {
            K = k; Lambda = lambda;
        }
        /// <summary>
        /// Gets or sets the value of the parameter k ∈ (0, +inf).
        /// </summary>
        public int K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter λ ∈ (0, +inf).
        /// </summary>
        public float Lambda
        {
            get
            {
                return this.lambda;
            }
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
                return this.k / this.lambda;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return this.k / Maths.Pow(this.lambda, 2.0f);
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                return new float[] { (this.k - 1f) / this.lambda };
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return 2.0f / Maths.Sqrt(k);
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
                return 6.0f / k;
            }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return (1 - k) * Special.DiGamma(k) + Maths.Log(Special.Gamma(k) / lambda) + k;
            }
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
                return 0;
            }
            return Maths.Pow(lambda, k) * Maths.Pow(x, k - 1) * Maths.Exp(-lambda * x) / (float)Special.Factorial(k - 1);
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Special.GammaIncomplete(k, lambda * x) / (float)Special.Factorial(k - 1);
        }
        #endregion
    }
}
