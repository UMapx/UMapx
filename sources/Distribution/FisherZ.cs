using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Fisher's Z-distribution.
    /// </summary>
    /// <remarks>
    /// The distribution is defined for positive degrees of freedom
    /// <c>d1 &gt; 0</c> and <c>d2 &gt; 0</c>. Its median has no closed form and is
    /// computed numerically. More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Fisher%27s_z-distribution
    /// </remarks>
    [Serializable]
    public class FisherZ : IDistribution
    {
        #region Private data
        private float d1;
        private float d2;
        #endregion

        #region FisherZ distribution
        /// <summary>
        /// Initializes the Fisher Z-distribution.
        /// </summary>
        /// <param name="d1">Degree of freedom d1 > 0</param>
        /// <param name="d2">Degree of freedom d2 > 0</param>
        public FisherZ(float d1, float d2)
        {
            D1 = d1; D2 = d2;
        }
        /// <summary>
        /// Gets or sets the degree of freedom d1 > 0.
        /// </summary>
        public float D1
        {
            get
            {
                return this.d1;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.d1 = value;
            }
        }
        /// <summary>
        /// Gets or sets the degree of freedom d2 > 0.
        /// </summary>
        public float D2
        {
            get
            {
                return this.d2;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.d2 = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                if (d2 <= 2f)
                    return float.NaN;

                return 0.5f * (Special.DiGamma(d1 / 2f) - Special.DiGamma(d2 / 2f) + Maths.Log(d2 / d1));
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (d2 <= 4f)
                    return float.NaN;

                return 0.25f * (Special.TriGamma(d1 / 2f) + Special.TriGamma(d2 / 2f));
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                return new float[] { 0f };
            }
        }
        /// <summary>
        /// Gets the median value (computed numerically).
        /// </summary>
        public float Median
        {
            get
            {
                IFloat f = (t) => Distribution(t) - 0.5f;

                float a = -10f, b = 10f;
                float fa = f(a), fb = f(b);

                while (fa * fb > 0f)
                {
                    a *= 2f; b *= 2f;
                    fa = f(a); fb = f(b);
                }

                for (int i = 0; i < 100; i++)
                {
                    float m = 0.5f * (a + b);
                    float fm = f(m);

                    if (Math.Abs(fm) < 1e-6f)
                        return m;

                    if (fa * fm < 0f)
                    {
                        b = m; fb = fm;
                    }
                    else
                    {
                        a = m; fa = fm;
                    }
                }

                return 0.5f * (a + b);
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
        /// Returns the value of differential entropy.
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
            get { return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity); }
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <remarks>
        /// The CDF is expressed through the regularized incomplete beta function
        /// after converting <paramref name="x"/> to the corresponding
        /// F-distribution variable.
        /// </remarks>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            float y = d1 * Maths.Exp(2f * x);
            if (float.IsPositiveInfinity(y)) return 1f;
            float u = y / (y + d2);
            return Special.BetaIncompleteRegularized(d1 * 0.5f, d2 * 0.5f, u);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <remarks>
        /// The density is given by:
        /// f(x; d1, d2) = 2^{-(d1 + d2)/2} · d1^{d1/2} · d2^{d2/2} · e^{d1·x} / (B(d1/2, d2/2) · (d1·e^{2x} + d2)^{(d1 + d2)/2}).
        /// </remarks>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            // helpers:

            float d12 = d1 / 2.0f;
            float d22 = d2 / 2.0f;
            float h = d12 + d22;

            // first equation:

            float a = Maths.Pow(d1, d12);
            float b = Maths.Pow(d2, d22);
            float c = 2f * a * b;
            float d = Special.Beta(d12, d22);
            float e = c / d;

            // second equation:

            float f = Maths.Exp(d1 * x);
            float g = d1 * Maths.Exp(2f * x) + d2;
            float j = f / Maths.Pow(g, h);

            // result of F(x, d1, d2):
            return e * j;
        }
        #endregion
    }
}
