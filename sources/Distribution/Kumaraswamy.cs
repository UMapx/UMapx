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
        private double a;
        private double b;
        #endregion

        #region Distribution
        /// <summary>
        /// Initializes the distribution of Kumarasva.
        /// </summary>
        /// <param name="a">Form parameter a > 0</param>
        /// <param name="b">Form parameter b > 0</param>
        public Kumaraswamy(double a, double b)
        {
            this.A = a;
            this.B = b;
        }
        /// <summary>
        /// Gets or sets form parameter a > 0.
        /// </summary>
        public double A
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
        public double B
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
        public double Mean
        {
            get
            {
                double num = b * Special.Gamma(1 + (1 / a)) * Special.Gamma(b);
                double den = Special.Gamma(1 + (1 / a) + b);

                return num / den;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                double alpha = momentGeneratingFunction(2, a, b);
                double beta = Math.Pow(momentGeneratingFunction(1, a, b), 2);
                return alpha - beta;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Pow(1 - Math.Pow(2, -1 / b), 1 / a);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {

                if ((a >= 1) && (b >= 1) && (a != 1 && b != 1))
                {
                    double num = a - 1;
                    double den = a * b - 1;
                    return Math.Pow(num / den, 1 / a);
                }

                return Double.NaN;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { return double.NaN; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(0, 1); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x > 1)
                return 1;

            if (x < 0)
                return 0;

            double xa = Math.Pow(x, a);
            return 1 - Math.Pow(1 - xa, b);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x > 1)
                return 0;

            if (x < 0)
                return 0;

            return a * b * Math.Pow(x, a - 1) * Math.Pow(1 - Math.Pow(x, a), b - 1);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        protected static double momentGeneratingFunction(int n, double a, double b)
        {
            return (b * Special.Beta(1.0 + ((double)n) / a, b));
        }
        #endregion
    }
}
