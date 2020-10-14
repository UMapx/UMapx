using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Gompertz distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gompertz_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Gompertz : IDistribution
    {
        #region Private data
        private double eta;
        private double b;
        #endregion

        #region Gompertz distribution
        /// <summary>
        ///Initializes the Gompertz distribution.
        /// </summary>
        /// <param name="eta">Form parameter η > 0</param>
        /// <param name="b">Scale parameter b > 0</param>
        public Gompertz(double eta, double b)
        {
            Eta = eta; B = b;
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter η > 0.
        /// </summary>
        public double Eta
        {
            get
            {
                return this.eta;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.eta = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter b > 0.
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
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
                if (eta >= 1)
                    return 0;

                return (1 / b) * Math.Log(1 / eta);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return (1.0 / b) * Math.Log((-1 / eta) * Math.Log(0.5) + 1);
            }
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
            get { return new RangeDouble(0, Double.PositiveInfinity); }
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

            double ebx = Math.Exp(b * x);
            return 1.0 - Math.Exp(-eta * (ebx - 1.0));
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

            double a1 = b * eta * Math.Exp(eta);
            double a2 = Math.Exp(b * x);
            double a3 = Math.Exp(-eta * Math.Exp(b * x));
            return a1 * a2 * a3;
        }
        #endregion
    }
}
