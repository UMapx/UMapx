using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the logistic distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Logistic_distribution
    /// </remarks>
    [Serializable]
    public class Logistic : IDistribution
    {
        #region Private data
        private float mu = 5;
        private float s = 2;
        #endregion

        #region Logistic components
        /// <summary>
        /// Initializes the logistic distribution.
        /// </summary>
        /// <param name="mu">Parameter μ</param>
        /// <param name="s">Parameter s (0, +inf]</param>
        public Logistic(float mu, float s)
        {
            Mu = mu; S = s;
        }
        /// <summary>
        /// Initializes the logistic distribution.
        /// </summary>
        public Logistic() { }
        /// <summary>
        /// Gets or sets the value of the parameter μ.
        /// </summary>
        public float Mu
        {
            get
            {
                return mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter s (0, +inf].
        /// </summary>
        public float S
        {
            get
            {
                return s;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.s = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return mu;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return (s * s * Maths.Pi * Maths.Pi) / 3.0f; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return mu;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get { return mu; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                return 1.2f;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get
            {
                return Maths.Log(s) + 2;
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            float z = (x - mu) / s;
            return 1.0f / (1 + Maths.Exp(-z));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            float z = (x - mu) / s;
            float num = Maths.Exp(-z);
            float a = (1 + num);
            float den = s * a * a;

            return num / den;
        }
        #endregion
    }
}
