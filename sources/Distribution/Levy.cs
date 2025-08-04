using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Levy distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/L%C3%A9vy_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Levy : IDistribution
    {
        #region Prviate data
        private float mu = 0;
        private float scale = 1;
        #endregion

        #region Levy components
        /// <summary>
        /// Initializes the Levy distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ</param>
        /// <param name="c">Scale factor (>0)</param>
        public Levy(float mu, float c)
        {
            Mu = mu; C = c;
        }
        /// <summary>
        /// Gets or sets the shift factor.
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
        /// Gets or sets the scale factor (> 0).
        /// </summary>
        public float C
        {
            get
            {
                return scale;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.scale = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(mu, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return float.PositiveInfinity; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return float.PositiveInfinity; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                if (mu == 0)
                {
                    return scale / 3.0f;
                }
                return float.NaN;
            }
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
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get
            {
                return (1.0f + 3.0f * Maths.Gamma + (float)Math.Log(16 * Math.PI * scale * scale)) / 2.0f;
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < mu)
            {
                return 0;
            }

            return Special.Erfc((float)Math.Sqrt(scale / (2 * (x - mu))));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < mu)
            {
                return 0;
            }
            float z = x - mu;
            float a = (float)Math.Sqrt(scale / (2.0 * Math.PI));
            float b = (float)Math.Exp(-(scale / (2 * z)));
            float c = (float)Math.Pow(z, 3.0 / 2.0);

            return a * b / c;
        }
        #endregion
    }
}
