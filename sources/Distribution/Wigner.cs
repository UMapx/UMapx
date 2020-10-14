using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Wiener semicircular distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Wigner : IDistribution
    {
        #region Private data
        private double r;
        #endregion;

        #region Wigner components
        /// <summary>
        /// Initializes the Wiener semicircular distribution.
        /// </summary>
        /// <param name="r">Radius</param>
        public Wigner(double r)
        {
            R = r;
        }
        /// <summary>
        /// Gets or sets the radius value.
        /// </summary>
        public double R
        {
            get
            {
                return this.r;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.r = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(-r, r);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return r * r / 4.0;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return 0;
            }
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
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return -1;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (Math.Abs(x) > r)
            {
                return double.NaN;
            }

            double r2 = r * r, x2 = x * x;
            double a = Math.Sqrt(r2 - x2);
            double b = 2 / (Math.PI * r2);
            return b * a;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (Math.Abs(x) > r)
            {
                return double.NaN;
            }

            double r2 = r * r, x2 = x * x;
            double a = Math.Sqrt(r2 - x2);
            double b = x / (Math.PI * r2);
            double c = Math.Asin(x / r) / Maths.Pi;
            return 0.5 + b * a + c;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(Maths.Pi * r) - 1.0 / 2;
            }
        }
        #endregion
    }
}
