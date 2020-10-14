using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Cauchy distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Cauchy_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Cauchy : IDistribution
    {
        #region Private data
        private double g = 0.5;
        private double x0 = 0;
        #endregion

        #region Caushi components
        /// <summary>
        /// Initializes the Cauchy distribution.
        /// </summary>
        public Cauchy() { }
        /// <summary>
        /// Initializes the Cauchy distribution.
        /// </summary>
        /// <param name="gamma">Scale factor (0, + inf)</param>
        /// <param name="x0">Shift coefficient</param>
        public Cauchy(double gamma, double x0)
        {
            Gamma = gamma;
            X0 = x0;
        }
        /// <summary>
        /// Gets or sets the value of the scale factor (0, + inf).
        /// </summary>
        public double Gamma
        {
            get
            {
                return this.g;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.g = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the shift coefficient.
        /// </summary>
        public double X0
        {
            get
            {
                return this.x0;
            }
            set
            {
                this.x0 = value;
            }
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
            get
            {
                return double.PositiveInfinity;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return x0;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return x0;
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
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return 1.0 / (Maths.Pi * g * (1.0 + Maths.Pow((x - x0) / g)));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            return 1.0 / Maths.Pi * Maths.Atg((x - x0) / g) + 0.5;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(4 * Maths.Pi * g);
            }
        }
        #endregion
    }
}
