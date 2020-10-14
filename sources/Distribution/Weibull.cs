using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Weibull distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Weibull_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Weibull : IDistribution
    {
        #region Private data
        private double l = 1;
        private double k = 1;
        #endregion

        #region Weibull components
        /// <summary>
        /// Initializes the Weibull distribution.
        /// </summary>
        public Weibull() { }
        /// <summary>
        /// Initializes the Weibull distribution.
        /// </summary>
        /// <param name="lambda">Scale factor (0, + inf)</param>
        /// <param name="k">Shape factor (0, + inf)</param>
        public Weibull(double lambda, double k)
        {
            Lambda = lambda;
        }
        /// <summary>
        /// Gets or sets the value of the scale factor (0, + inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.l = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the form factor (0, + inf).
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                this.k = Maths.Max(0.0000001, value);
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
                return l * Special.Gamma(1.0 + 1.0 / k);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return l * l * (Special.Gamma(1.0 + 2.0 / k) - Special.Gamma(1.0 + 1.0 / k));
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return l * Maths.Pow(k - 1, 1.0 / k) / Maths.Pow(k, 1.0 / k);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return l * Maths.Pow(Maths.Log(2), 1.0 / k);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (Special.Gamma(1.0 + 3.0 / k) * Maths.Pow(l, 3) - 3 * Mean * Special.Gamma(1 + 2.0 / k) * Maths.Pow(l, 2) + 2 * Maths.Pow(Mean, 3)) / (Variance * Maths.Sqrt(Variance));
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return (Special.Gamma(1.0 + 4.0 / k) * Maths.Pow(l, 4) - 4 * Mean * Special.Gamma(1 + 3.0 / k) * Maths.Pow(l, 3) + 6 * Maths.Pow(Mean, 2) * Math.Pow(l, 2) * Special.Gamma(1.0 + 2.0 / k) - 3 * Maths.Pow(Mean, 4)) / (Variance * Variance);
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
            return (k / l) * Maths.Pow(x / l, k - 1) * Maths.Exp(-Maths.Pow(x / l, k));
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
            return 1 - Maths.Exp(-Maths.Pow(x / l, k));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Maths.Gamma * (1.0 - 1.0 / k) + Maths.Pow(l / k, k) + Maths.Log(l / k);
            }
        }
        #endregion
    }
}
