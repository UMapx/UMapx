﻿using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Gamma-distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gamma_distribution
    /// </remarks>
    [Serializable]
    public class Gamma : IDistribution
    {
        #region Private data
        private float thetta = 1;
        private float k = 1;
        #endregion

        #region Gamma components
        /// <summary>
        /// Initializes the Gamma-distribution.
        /// </summary>
        public Gamma() { }
        /// <summary>
        /// Initializes the Gamma-distribution.
        /// </summary>
        /// <param name="thetta">Parameter θ (0, +inf)</param>
        /// <param name="k">Parameter k (0, +inf)</param>
        public Gamma(float thetta, float k)
        {
            Thetta = thetta; K = k;
        }
        /// <summary>
        /// Gets or sets the parameter θ (0, +inf).
        /// </summary>
        public float Thetta
        {
            get
            {
                return this.thetta;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.thetta = value;
            }
        }
        /// <summary>
        /// Gets or sets the parameter k (0, +inf).
        /// </summary>
        public float K
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
                return thetta * k;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                return k * thetta * thetta;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                if (k >= 1)
                {
                    return (k - 1) * thetta;
                }
                return 0;
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
                return 2 / Maths.Sqrt(k);
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
            return Maths.Pow(x, k - 1) * Maths.Exp(-x / thetta) / (Special.Gamma(k) * Maths.Pow(thetta, k));
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
            return Special.GammaP(k, x / thetta);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                return k + Maths.Log(thetta) + Special.LogGamma(k) + (1f - k) * Special.DiGamma(k);
            }
        }
        #endregion
    }
}
