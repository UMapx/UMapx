using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Levy distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/L%C3%A9vy_distribution
    /// </remarks>
    [Serializable]
    public class Levy : IDistribution
    {
        #region Private data
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
                    throw new ArgumentException("Invalid argument value");

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
                return mu + scale / 3.0f;
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
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the value of differential entropy (h = 1 + ln√(2πc)).
        /// </summary>
        public float Entropy
        {
            get
            {
                return 1f + 0.5f * Maths.Log(2f * Maths.Pi * scale);
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

            return Special.Erfc(Maths.Sqrt(scale / (2 * (x - mu))));
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
            float a = Maths.Sqrt(scale / (2.0f * Maths.Pi));
            float b = Maths.Exp(-(scale / (2 * z)));
            float c = Maths.Pow(z, 3.0f / 2.0f);

            return a * b / c;
        }
        #endregion
    }
}
