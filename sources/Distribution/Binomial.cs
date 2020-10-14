using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the binomial distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Binomial_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Binomial : IDistribution
    {
        #region Private data
        private double n = 20;
        private double p = 0.5;
        private double q = 0.5;
        #endregion

        #region Binomial components
        /// <summary>
        /// Initializes the binomial distribution.
        /// </summary>
        public Binomial() { }
        /// <summary>
        /// Initializes the binomial distribution.
        /// </summary>
        /// <param name="n">Number of experiments (>0)</param>
        /// <param name="p">Probability of success [0, 1]</param>
        public Binomial(double n, double p)
        {
            N = n; P = p;
        }
        /// <summary>
        /// Gets or sets number of experiments.
        /// </summary>
        public double N
        {
            get
            {
                return n;
            }
            set
            {
                this.n = value;
            }
        }
        /// <summary>
        /// Gets or sets probability of success [0, 1].
        /// </summary>
        public double P
        {
            get
            {
                return p;
            }
            set
            {
                if (value > 1 || value < 0)
                    throw new Exception("Invalid argument value");

                this.p = value;
                this.q = 1.0 - p;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, this.n);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return n * p;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return n * p * q;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                double test = (n + 1) * p;

                if (test <= 0 || (int)test != test)
                    return Math.Floor(test);

                if (test <= n)
                    return test;

                if (test == n + 1)
                    return n;

                return Double.NaN;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Floor(n * p);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (q - p) / Math.Sqrt(n * p * q);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return (1 - 6 * p * q) / Variance;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0 || x > n)
            {
                return 0;
            }

            double a = Special.LogBinomial(n, x);
            double b = x == 0 ? 0 : x * Math.Log(p);
            double c = (n - x);
            double d = Math.Log(1 - p);
            double log = a + b + c * d;

            return Math.Exp(log);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
                return 0;
            if (x >= n)
                return 1;

            double a = n - x;
            double b = x + 1;
            return Special.Beta(a, b, q);
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
