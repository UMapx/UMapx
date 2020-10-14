using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Levy distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/L%C3%A9vy_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Levy : IDistribution
    {
        #region Prviate data
        private double mu = 0;
        private double scale = 1;
        #endregion

        #region Levy components
        /// <summary>
        /// Initializes the Levy distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ</param>
        /// <param name="c">Scale factor (>0)</param>
        public Levy(double mu, double c)
        {
            Mu = mu; C = c;
        }
        /// <summary>
        /// Gets or sets the shift factor.
        /// </summary>
        public double Mu
        {
            get
            {
                return mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor (> 0).
        /// </summary>
        public double C
        {
            get
            {
                return scale;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.scale = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(mu, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return Double.PositiveInfinity; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return Double.PositiveInfinity; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (mu == 0)
                {
                    return scale / 3.0;
                }
                return Double.NaN;
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
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                return (1.0 + 3.0 * Maths.Gamma + Math.Log(16 * Math.PI * scale * scale)) / 2.0;
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < mu)
            {
                return 0;
            }

            return Special.Erfc(Math.Sqrt(scale / (2 * (x - mu))));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < mu)
            {
                return 0;
            }
            double z = x - mu;
            double a = Math.Sqrt(scale / (2.0 * Math.PI));
            double b = Math.Exp(-(scale / (2 * z)));
            double c = Math.Pow(z, 3.0 / 2.0);

            return a * b / c;
        }
        #endregion
    }
}
