using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the triangular distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Triangular_distribution
    /// </remarks>
    [Serializable]
    public class Triangular : IDistribution
    {
        #region Private data
        private float a;
        private float b;
        private float c;
        #endregion

        #region Triangular components
        /// <summary>
        /// Initializes the triangular distribution.
        /// </summary>
        public Triangular() { }
        /// <summary>
        /// Initializes the triangular distribution.
        /// </summary>
        /// <param name="a">Parameter a ∈ (-inf, +inf)</param>
        /// <param name="b">Parameter b ∈ (-inf, +inf)</param>
        /// <param name="c">Parameter c ∈ (-inf, +inf)</param>
        public Triangular(float a, float b, float c)
        {
            A = a; B = b; C = c;
        }
        /// <summary>
        /// Gets or sets the value of the parameter a ∈ (-inf, +inf).
        /// </summary>
        public float A
        {
            get
            {
                return a;
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter b ∈ (-inf, +inf).
        /// </summary>
        public float B
        {
            get
            {
                return b;
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter c ∈ (-inf, +inf).
        /// </summary>
        public float C
        {
            get
            {
                return c;
            }
            set
            {
                this.c = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(a, b);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { return (a + b + c) / 3.0f; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { return (a * a + b * b + c * c - a * b - a * c - b * c) / 18; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        /// <remarks>
        /// The <c>c ≥ (a + b) / 2</c> branch returns the median on the rising side of the PDF;
        /// otherwise the median lies on the falling side.
        /// </remarks>
        public float Median
        {
            get
            {
                if (c >= (a + b) / 2.0f)
                {
                    return a + Maths.Sqrt((b - a) * (c - a) / 2.0f);
                }

                return b - Maths.Sqrt((b - a) * (b - c) / 2.0f);
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get { return new float[] { c }; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                float k1 = (a + b - 2 * c) * (2 * a - b - c) * (a - 2 * b + c);
                float k2 = 5 * (a * a + b * b + c * c - a * b - a * c - b * c);
                return Maths.Sqrt2 * k1 / Maths.Pow(k2, 1.5f);
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
                return -0.6f;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get
            {
                return 0.5f + Maths.Log((b - a) / 2);
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < a)
                return 0;

            if (x >= a && x <= c)
                return ((x - a) * (x - a)) / ((b - a) * (c - a));

            if (x > c && x <= b)
                return 1 - ((b - x) * (b - x)) / ((b - a) * (b - c));

            return 1;
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

            if (x >= a && x <= c)
                return (2 * (x - a)) / ((b - a) * (c - a));

            if (x > c && x <= b)
                return (2 * (b - x)) / ((b - a) * (b - c));

            return 0;
        }
        #endregion
    }
}
