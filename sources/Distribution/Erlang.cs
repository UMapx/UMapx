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
        private double lambda = 0.5;
        #endregion

        #region Erlang distribution
        /// <summary>
        /// Initializes the distribution of Erlang.
        /// </summary>
        /// <param name="k">Form parameter k ∈ (0, +inf)</param>
        /// <param name="lambda">λ-parameter λ ∈ (0, +inf)</param>
        public Erlang(int k, double lambda)
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
        public double Lambda
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
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return this.k / this.lambda;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return this.k / Math.Pow(this.lambda, 2.0);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return (this.k - 1) / this.lambda;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2.0 / Math.Sqrt(k);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6.0 / k;
            }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return (1 - k) * Special.DiGamma(k) + Math.Log(Special.Gamma(k) / lambda) + k;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return double.NaN;
            }
            return Math.Pow(lambda, k) * Math.Pow(x, k - 1) * Math.Exp(-lambda * x) / Special.Factorial(k - 1);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return double.NaN;
            }
            return Special.GammaIncomplete(k, lambda * x) / Special.Factorial(k - 1);
        }
        #endregion
    }
}
