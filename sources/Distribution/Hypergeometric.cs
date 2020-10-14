using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the hypergeometric distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hypergeometric_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Hypergeometric : IDistribution
    {
        #region Private data
        private double n = 30;
        private double k = 20;
        private double d = 20;
        #endregion

        #region Hypergeometric components
        /// <summary>
        /// Initializes the hypergeometric distribution.
        /// </summary>
        public Hypergeometric() { }
        /// <summary>
        /// Initializes the hypergeometric distribution.
        /// </summary>
        /// <param name="n">Parameter N [0, +inf]</param>
        /// <param name="k">Parameter D [0, N]</param>
        /// <param name="d">Parameter K [0, N]</param>
        public Hypergeometric(double n, double k, double d)
        {
            N = n; K = k; D = d;
        }
        /// <summary>
        /// Gets or sets the value of the parameter N [0, +inf].
        /// </summary>
        public double N
        {
            get
            {
                return n;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.n = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter D [0, N].
        /// </summary>
        public double D
        {
            get
            {
                return d;
            }
            set
            {
                if (value < 0 || value > N)
                    throw new Exception("Invalid argument value");

                this.d = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter k [0, N].
        /// </summary>
        public double K
        {
            get
            {
                return k;
            }
            set
            {
                if (value < 0 || value > N)
                    throw new Exception("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, k);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return k * d / n;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return k * (d / n) * ((n - d) / n) * ((n - k) / (n - 1.0)); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                double num = (k + 1) * (d + 1);
                double den = n + 2;
                return Math.Floor(num / den);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                double k1 = (n - 2 * d) * Math.Pow(n - 1, 0.5) * (n - 2 * k);
                double k2 = (n * d * (n - d) * Math.Pow(n - k, 0.5) * (n - 2));
                return k1 / k2;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                double n2 = n * n;
                double k1 = (n2 * (n - 1)) / (k * (n - 2) * (n - 3) * (n - k));
                double k2 = (n * (n + 1) - 6 * n * (n - k)) / (d * (n - d));
                double k3 = (3 * n * (n - k) * (n + 6)) / n2 - 6;
                return k1 * (k2 + k3);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < Math.Max(0, k + d - n) || x > Math.Min(d, k))
            {
                return 0;
            }

            double a = Special.Binomial(d, x);
            double b = Special.Binomial(n - d, k - x);
            double c = Special.Binomial(n, k);
            return (a * b) / c;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            double sum = 0;
            int k = (int)x;
            for (int i = 0; i <= k; i++)
            {
                sum += Function(i);
            }

            return sum;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
}
