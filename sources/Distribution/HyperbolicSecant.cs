﻿using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the hyperbolic secant distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class HyperbolicSecant : IDistribution
    {
        #region Distribution
        /// <summary>
        /// Initializes the hyperbolic secant distribution.
        /// </summary>
        public HyperbolicSecant() { }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return 1; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get { return 2; }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity); }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        public float Entropy
        {
            get { return (4.0f / MathsF.Pi) * MathsF.G; }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            float angle = (float)Math.Atan(Math.Exp(x * Math.PI / 2.0));
            return 2 * angle / MathsF.Pi;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            return 0.5f * MathsF.Sech(x * (MathsF.Pi / 2.0f));
        }
        #endregion
    }
}
