using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Poisson distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Poisson_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Poisson : IDistribution
    {
        #region Private data
        private double l = 1;
        #endregion

        #region Poisson components
        /// <summary>
        /// Initializes the Poisson distribution.
        /// </summary>
        public Poisson() { }
        /// <summary>
        /// Initializes the Poisson distribution.
        /// </summary>
        /// <param name="lambda">Parameter λ (0, +inf)</param>
        public Poisson(double lambda)
        {
            Lambda = lambda;
        }
        /// <summary>
        /// Gets or sets the value of the parameter λ (0, +inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.l = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return l;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return l;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return Math.Floor(l);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Floor(l + 0.333 - 0.02 / l);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return Math.Pow(l, -0.5);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return Math.Pow(l, -1.0);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Math.Exp(-l) * Math.Pow(l, x) / Special.Factorial(x);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Special.GammaQ(x + 1, l);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return l * (1 - Math.Log(l)) + Math.Exp(-l) * Row(l);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="l"></param>
        /// <returns></returns>
        private double Row(double l)
        {
            double sum = 0;
            int k, n = 20;
            double fac;

            for (k = 0; k < n; k++)
            {
                fac = Special.Factorial(k);
                sum += Math.Pow(l, k) * Math.Log(fac) / fac;
            }

            return sum;
        }
        #endregion
    }
}
