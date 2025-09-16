using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the trapezoidal distribution.
    /// </summary>
    [Serializable]
    public class Trapezoidal : IDistribution
    {
        #region Private data
        private float a;
        private float b;
        private float c;
        private float d;
        private float n1;
        private float n3;
        private float alpha;
        private float constant;
        #endregion

        #region Constructor
        /// <summary>
        /// Initializes the trapezoidal distribution with default slopes and boundary ratio.
        /// </summary>
        /// <param name="a">Minimum value.</param>
        /// <param name="b">Beginning of the stability region.</param>
        /// <param name="c">End of the stability region.</param>
        /// <param name="d">Maximum value.</param>
        public Trapezoidal(float a, float b, float c, float d)
            : this(a, b, c, d, 2f, 2f, 1f)
        {
        }
        /// <summary>
        /// Initializes the trapezoidal distribution with specified slopes.
        /// </summary>
        /// <param name="a">Minimum value.</param>
        /// <param name="b">Beginning of the stability region.</param>
        /// <param name="c">End of the stability region.</param>
        /// <param name="d">Maximum value.</param>
        /// <param name="n1">Growth slope.</param>
        /// <param name="n3">Decay slope.</param>
        public Trapezoidal(float a, float b, float c, float d, float n1, float n3)
            : this(a, b, c, d, n1, n3, 1f)
        {
        }
        /// <summary>
        /// Initializes the trapezoidal distribution with the specified parameters.
        /// </summary>
        /// <param name="a">Minimum value.</param>
        /// <param name="b">Beginning of the stability region.</param>
        /// <param name="c">End of the stability region.</param>
        /// <param name="d">Maximum value.</param>
        /// <param name="n1">Growth slope.</param>
        /// <param name="n3">Decay slope.</param>
        /// <param name="alpha">Boundary ratio.</param>
        public Trapezoidal(float a, float b, float c, float d, float n1, float n3, float alpha)
        {
            if (a > b)
            {
                throw new ArgumentOutOfRangeException(nameof(b), "Parameter b must be greater than or equal to a");
            }
            if (b > c)
            {
                throw new ArgumentOutOfRangeException(nameof(c), "Parameter c must be greater than or equal to b");
            }
            if (d < c)
            {
                throw new ArgumentOutOfRangeException(nameof(d), "Parameter d must be greater than or equal to c");
            }
            if (d <= a)
            {
                throw new ArgumentOutOfRangeException(nameof(d), "Parameter d must be greater than a");
            }
            if (n1 <= 0)
            {
                throw new ArgumentOutOfRangeException(nameof(n1), "Slope n1 must be positive");
            }
            if (n3 <= 0)
            {
                throw new ArgumentOutOfRangeException(nameof(n3), "Slope n3 must be positive");
            }
            if (alpha <= 0)
            {
                throw new ArgumentOutOfRangeException(nameof(alpha), "Alpha must be positive");
            }

            this.a = a;
            this.b = b;
            this.c = c;
            this.d = d;
            this.n1 = n1;
            this.n3 = n3;
            this.alpha = alpha;

            double num = 2.0 * n1 * n3;
            double den = 2.0 * alpha * (b - a) * n3 + (alpha + 1.0) * (c - b) * n1 * n3 + 2.0 * (d - c) * n1;
            this.constant = (float)(num / den);
        }
        #endregion

        #region Distribution properties
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support => new RangeFloat(a, d);
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                double expectationX1 = (a + n1 * b) / (n1 + 1.0);
                double expectationX3 = (n3 * c + d) / (n3 + 1.0);

                double num = (-2.0 / 3.0) * (alpha - 1.0) * (Math.Pow(c, 3) - Math.Pow(b, 3))
                    + (alpha * c - b) * (Math.Pow(c, 2) - Math.Pow(b, 2));
                double den = Math.Pow(c - b, 2) * (alpha + 1.0);
                double expectationX2 = num / den;

                num = (2.0 * alpha * (b - a) * n3 * expectationX1)
                    + (n1 * n3 * expectationX2)
                    + (2.0 * (d - c) * n1 * expectationX3);

                den = (2.0 * alpha * (b - a) * n3)
                    + ((alpha + 1.0) * (c - b) * n1 * n3)
                    + (2.0 * (d - c) * n1);

                return (float)(num / den);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                double expectationX1_2 = (2.0 * a * a + 2.0 * n1 * a * b + n1 * (n1 + 1.0) * b * b)
                    / ((n1 + 2.0) * (n1 + 1.0));

                double num = -0.5 * (alpha - 1.0) * (Math.Pow(c, 4) - Math.Pow(b, 4))
                    + (2.0 / 3.0) * (alpha * c - b) * (Math.Pow(c, 3) - Math.Pow(b, 3));
                double den = Math.Pow(c - b, 2) * (alpha + 1.0);
                double expectationX2_2 = num / den;

                double expectationX3_2 = (2.0 * d * d + 2.0 * n3 * c * d + n3 * (n3 + 1.0) * c * c)
                    / ((n3 + 2.0) * (n3 + 1.0));

                double common = 2.0 * alpha * (b - a) * n3 + (alpha + 1.0) * (c - b) * n1 * n3 + 2.0 * (d - c) * n1;

                double x = (2.0 * alpha * (b - a) * n3) / common;
                double y = (n1 * n3) / common;
                double z = (2.0 * (d - c) * n1) / common;

                return (float)(x * expectationX1_2 + y * expectationX2_2 + z * expectationX3_2);
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
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get { return new float[] { float.NaN }; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness => float.NaN;
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        public float Excess => float.NaN;
        /// <summary>
        /// Gets the value of differential entropy.
        /// </summary>
        public float Entropy => float.NaN;
        #endregion

        #region Methods
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        public float Function(float x)
        {
            if (x < a || x > d)
            {
                return 0f;
            }

            if (x < b)
            {
                double ratio = (x - a) / (b - a);
                return (float)(constant * alpha * Math.Pow(ratio, n1 - 1.0));
            }

            if (x < c)
            {
                double term = ((alpha - 1.0) * (c - x) / (c - b)) + 1.0;
                return (float)(constant * term);
            }

            if (x < d)
            {
                double ratio = (d - x) / (d - c);
                return (float)(constant * Math.Pow(ratio, n3 - 1.0));
            }

            if (Math.Abs(x - d) < 1e-6f)
            {
                return 0f;
            }

            return 0f;
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        public float Distribution(float x)
        {
            if (x <= a)
            {
                return 0f;
            }

            if (x < b)
            {
                double num = 2.0 * alpha * (b - a) * n3;
                double den = 2.0 * alpha * (b - a) * n3 + (alpha + 1.0) * (c - b) * n1 * n3 + 2.0 * (d - c) * n1;
                double p = Math.Pow((x - a) / (b - a), n1);
                return (float)((num / den) * p);
            }

            if (x < c)
            {
                double num = 2.0 * alpha * (b - a) * n3
                    + 2.0 * (x - b) * n1 * n3 * (1.0 + ((alpha - 1.0) / 2.0) * ((2.0 * c - b - x) / (c - b)));
                double den = 2.0 * alpha * (b - a) * n3 + (alpha + 1.0) * (c - b) * n1 * n3 + 2.0 * (d - c) * n1;
                return (float)(num / den);
            }

            if (x < d)
            {
                double num = 2.0 * (d - c) * n1;
                double den = 2.0 * alpha * (b - a) * n3 + (alpha + 1.0) * (c - b) * n1 * n3 + 2.0 * (d - c) * n1;
                double p = Math.Pow((d - x) / (d - c), n3);
                return (float)(1.0 - (num / den) * p);
            }

            return 1f;
        }
        #endregion
    }
}
