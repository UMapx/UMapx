using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the U-quadratic distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/U-quadratic_distribution
    /// </remarks>
    [Serializable]
    public class UQuadratic : IDistribution
    {
        #region Private data
        float a;
        float b;
        float alpha;
        float beta;

        /// <summary>
        /// Updates coefficients when parameters change.
        /// </summary>
        private void UpdateCoefficients()
        {
            this.alpha = 12 / Maths.Pow(b - a, 3);
            this.beta = (b + a) / 2;
        }
        #endregion

        #region UQuadratic components
        /// <summary>
        /// Initializes the U-quadratic distribution.
        /// </summary>
        /// <param name="a">Parameter a ∈ (0, +inf)</param>
        /// <param name="b">Parameter b ∈ (a, +inf)</param>
        public UQuadratic(float a, float b)
        {
            A = a; B = b;
            UpdateCoefficients();
        }
        /// <summary>
        /// Gets or sets the value of the parameter a ∈ (0, +inf).
        /// </summary>
        public float A
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.a = value;
                UpdateCoefficients();
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter b ∈ (a, +inf).
        /// </summary>
        public float B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (value <= a)
                    throw new ArgumentException("The value of parameter b must be greater than a");

                this.b = value;
                UpdateCoefficients();
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                return new float[] { a, b };
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return (a + b) / 2.0f; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { return (a + b) / 2.0f; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return (3f / 20f) * Maths.Pow(b - a, 2f); }
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
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get
            {
                return (3f / 112f * Maths.Pow(b - a, 4f)) / Maths.Pow(Variance, 2f) - 3f;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(a, b); }
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < a)
                return 0;

            if (x > b)
                return 1;

            return (alpha / 3) * (float)(Math.Pow(x - beta, 3) + Maths.Pow(beta - a, 3));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < a)
                return 0;

            if (x > b)
                return 0;

            return alpha * Maths.Pow(x - beta, 2);
        }
        #endregion
    }
}
