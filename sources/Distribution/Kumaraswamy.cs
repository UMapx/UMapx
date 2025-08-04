using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the distribution of Kumaraswa.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Kumaraswamy_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Kumaraswamy : IDistribution
    {
        #region Private data
        private float a;
        private float b;
        #endregion

        #region Distribution
        /// <summary>
        /// Initializes the distribution of Kumarasva.
        /// </summary>
        /// <param name="a">Form parameter a > 0</param>
        /// <param name="b">Form parameter b > 0</param>
        public Kumaraswamy(float a, float b)
        {
            this.A = a;
            this.B = b;
        }
        /// <summary>
        /// Gets or sets form parameter a > 0.
        /// </summary>
        public float A
        {
            get
            {
                return a;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets form parameter b > 0.
        /// </summary>
        public float B
        {
            get
            {
                return b;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.b = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                float num = b * Special.Gamma(1 + (1 / a)) * Special.Gamma(b);
                float den = Special.Gamma(1 + (1 / a) + b);

                return num / den;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                float alpha = momentGeneratingFunction(2, a, b);
                float beta = (float)Math.Pow(momentGeneratingFunction(1, a, b), 2);
                return alpha - beta;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return (float)Math.Pow(1 - Math.Pow(2, -1 / b), 1 / a);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {

                if ((a >= 1) && (b >= 1) && (a != 1 && b != 1))
                {
                    float num = a - 1;
                    float den = a * b - 1;
                    return (float)Math.Pow(num / den, 1 / a);
                }

                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get { return float.NaN; }
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
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(0, 1); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x > 1)
                return 1;

            if (x < 0)
                return 0;

            float xa = (float)Math.Pow(x, a);
            return 1 - (float)Math.Pow(1 - xa, b);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x > 1)
                return 0;

            if (x < 0)
                return 0;

            return a * b * (float)Math.Pow(x, a - 1) * (float)Math.Pow(1 - Math.Pow(x, a), b - 1);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        protected static float momentGeneratingFunction(int n, float a, float b)
        {
            return (b * Special.Beta(1.0f + ((float)n) / a, b));
        }
        #endregion
    }
}
