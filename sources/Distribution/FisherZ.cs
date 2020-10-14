using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Fisher's Z-distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Fisher%27s_z-distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FisherZ : IDistribution
    {
        #region Private data
        private double d1;
        private double d2;
        #endregion

        #region FisherZ distribution
        /// <summary>
        /// Initializes the Fisher Z-distribution.
        /// </summary>
        /// <param name="d1">Degree of freedom d1 > 0</param>
        /// <param name="d2">Degree of freedom d2 > 0</param>
        public FisherZ(double d1, double d2)
        {
            D1 = d1; D2 = d2;
        }
        /// <summary>
        /// Gets or sets the degree of freedom d1 > 0.
        /// </summary>
        public double D1
        {
            get
            {
                return this.d1;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.d1 = value;
            }
        }
        /// <summary>
        /// Gets or sets the degree of freedom d2 > 0.
        /// </summary>
        public double D2
        {
            get
            {
                return this.d2;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.d2 = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
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
        /// Returns the value of differential entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(Double.NegativeInfinity, Double.PositiveInfinity); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            // helpers:

            double d12 = d1 / 2.0;
            double d22 = d2 / 2.0;

            // first equation:

            double a = Math.Pow(d1, d12);
            double b = Math.Pow(d2, d22);
            double c = 2 * a * b;
            double d = Special.Beta(d12, d22);
            double e = c / d;

            // second equation:

            double f = Math.Exp(d1 * x);
            double g = d1 * Math.Exp(2 * x) + d2;
            double h = d12 + d22;
            double j = f / Math.Pow(g, h);

            // result of F(x, d1, d2):
            return e * j;
        }
        #endregion
    }
}
