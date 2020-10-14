using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the distribution of Choi Williams.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Choi%E2%80%93Williams_distribution_function
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ChoiWilliams : IDistribution
    {
        #region Private data
        private double a;
        #endregion

        #region Choi-Williams components
        /// <summary>
        /// Initializes the Choi-Williams distribution.
        /// </summary>
        /// <param name="a">Coefficient</param>
        public ChoiWilliams(double a = 0.001)
        {
            A = a;
        }
        /// <summary>
        /// Gets or sets the coefficient value.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
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
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
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
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the kernel density function.
        /// </summary>
        /// <param name="eta">Argument</param>
        /// <param name="tau">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double eta, double tau)
        {
            double ksi = eta * tau;
            return Math.Exp(-a * ksi * ksi);
        }
        /// <summary>
        /// Returns the value of the kernel distribution function.
        /// </summary>
        /// <param name="t">Argument</param>
        /// <param name="tau">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double t, double tau)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
}
