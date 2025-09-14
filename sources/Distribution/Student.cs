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
                if (value <= 0)
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
        /// <remarks>
        /// For degrees of freedom ≤ 1, the mean is undefined and returns <see cref="float.NaN"/>.
        /// </remarks>
        public float Mean
        {
            get => degrees > 1 ? 0f : float.NaN;
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
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// For degrees of freedom in the interval (2, 4], the excess kurtosis is infinite and
        /// returns <see cref="float.PositiveInfinity"/>. For degrees of freedom ≤ 2, the
        /// excess kurtosis is undefined and returns <see cref="float.NaN"/>.
        /// </remarks>
        public float Excess
        {
            get
            {
                if (degrees > 4)
                    return 6f / (degrees - 4);
                if (degrees > 2)
                    return float.PositiveInfinity;
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
            float u = Special.BetaIncompleteRegularized(v / 2.0f, 0.5f, v / (x * x + v));
            return x >= 0 ? 1.0f - 0.5f * u : 0.5f * u;
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
            float num = Special.LogGamma((this.degrees + 1) / 2.0f);
            float den = 0.5f * Maths.Log(this.degrees * Maths.Pi) + Special.LogGamma(this.degrees / 2.0f);
            float lambda = num - den;
            return lambda - (degrees + 1) / 2.0f * Maths.Log(1 + (x * x) / degrees);
        }
        #endregion
    }
}
