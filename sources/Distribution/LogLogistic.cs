using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the log-logistic distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Log-logistic_distribution
    /// </remarks>
    [Serializable]
    public class LogLogistic : IDistribution
    {
        #region Private data
        private float a = 1;
        private float b = 1;
        #endregion

        #region LogLogistic components
        /// <summary>
        /// Initializes the log-logistic distribution.
        /// </summary>
        public LogLogistic() { }
        /// <summary>
        /// Initializes the log-logistic distribution.
        /// </summary>
        /// <param name="a">Parameter a</param>
        /// <param name="b">Parameter b</param>
        public LogLogistic(float a, float b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the value of parameter a.
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
                    throw new ArgumentException("Invalid argument value");

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of parameter b.
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
                    throw new ArgumentException("Invalid argument value");

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
                return new RangeFloat(0, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                if (b > 1)
                {
                    return a * (float)Math.Pow((b - 1) / (b + 1), 1 / b);
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
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= 0)
            {
                return 0;
            }

            return (b / a) * (float)Math.Pow(x / a, b - 1) / (float)Math.Pow(1.0f + Math.Pow(x / a, b), 2);
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
                return 0;
            }

            return 1.0f / (1 + (float)Math.Pow(x / a, -b));
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
