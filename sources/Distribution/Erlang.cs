using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the distribution of Erlang.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Erlang_distribution
    /// </remarks>
    /// </summary>
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
                    throw new Exception("Invalid argument value");

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
                    throw new Exception("Invalid argument value");

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
                return this.k / (float)Math.Pow(this.lambda, 2.0);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return (this.k - 1) / this.lambda;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return 2.0f / (float)Math.Sqrt(k);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
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
                return (1 - k) * Special.DiGamma(k) + (float)Math.Log(Special.Gamma(k) / lambda) + k;
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
                return float.NaN;
            }
            return (float)Math.Pow(lambda, k) * (float)Math.Pow(x, k - 1) * (float)Math.Exp(-lambda * x) / Special.Factorial(k - 1);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < 0)
            {
                return float.NaN;
            }
            return Special.GammaIncomplete(k, lambda * x) / Special.Factorial(k - 1);
        }
        #endregion
    }
}
