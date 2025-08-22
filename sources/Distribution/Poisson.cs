using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Poisson distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Poisson_distribution
    /// </remarks>
    /// </summary>
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
        /// <param name="lambda">Parameter λ (0, +inf)</param>
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
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

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
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return (float)Math.Floor(l);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return (float)Math.Floor(l + 0.333 - 0.02 / l);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return (float)Math.Pow(l, -0.5);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                return (float)Math.Pow(l, -1.0);
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
            return (float)Math.Exp(-l) * (float)Math.Pow(l, x) / (float)Special.Factorial(x);
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
                return 0;
            }
            return Special.GammaQ(x + 1, l);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return l * (1 - (float)Math.Log(l)) + (float)Math.Exp(-l) * Row(l);
            }
        }
        /// <summary>
        /// Calculate row.
        /// </summary>
        /// <param name="l">Value</param>
        /// <returns>Value</returns>
        private float Row(float l)
        {
            float sum = 0;
            int k, n = 20;
            float fac;

            for (k = 0; k < n; k++)
            {
                fac = (float)Special.Factorial(k);
                sum += (float)Math.Pow(l, k) * (float)Math.Log(fac) / fac;
            }

            return sum;
        }
        #endregion
    }
}
