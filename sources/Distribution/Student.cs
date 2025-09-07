using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Student's distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Student%27s_t-distribution
    /// </remarks>
    [Serializable]
    public class Student : IDistribution
    {
        #region Private data
        private float lambda;
        private float degrees;
        #endregion

        #region Student components
        /// <summary>
        /// Initializes the Student's distribution.
        /// </summary>
        /// <param name="n">Degrees of freedom n ∈ (0, +inf)</param>
        public Student(float n)
        {
            this.N = n;
            float num = Special.LogGamma((n + 1) / 2.0f);
            float den = 0.5f * Maths.Log(n * Maths.Pi) + Special.LogGamma(n / 2.0f);
            this.lambda = num - den;
        }
        /// <summary>
        /// Gets or sets degrees of freedom n ∈ (0, +inf).
        /// </summary>
        public float N
        {
            get
            {
                return this.degrees;
            }
            set
            {
                if (value < 1)
                    throw new ArgumentException("Invalid argument value");

                this.degrees = value;
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
            get { return 0; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (degrees > 2)
                    return degrees / (degrees - 2);
                else if (degrees > 1)
                    return float.PositiveInfinity;
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                if (degrees > 3)
                {
                    return 0;
                }
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get
            {
                if (degrees > 4)
                {
                    return 6.0f / (degrees - 4);
                }
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get
            {
                float a = Special.DiGamma((1 + degrees) / 2.0f);
                float b = Special.DiGamma(degrees / 2.0f);
                float c = (degrees + 1) / 2.0f * (a - b);
                float d = Maths.Sqrt(degrees) * Special.Beta(degrees / 2.0f, 0.5f);

                return c + Maths.Log(d);
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            float v = degrees;
            float sqrt = Maths.Sqrt(x * x + v);
            float u = (x + sqrt) / (2 * sqrt);
            return Special.BetaIncompleteRegularized(v / 2.0f, v / 2.0f, u);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            return Maths.Exp(LogFunction(x));
        }
        /// <summary>
        /// Computes the natural logarithm of the probability density function.
        /// </summary>
        /// <param name="x">Input value</param>
        /// <returns>Logarithm of the density</returns>
        private float LogFunction(float x)
        {
            return lambda - ((degrees + 1) / 2.0f) * Maths.Log(1 + (x * x) / degrees);
        }
        #endregion
    }
}
