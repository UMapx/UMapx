using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Fisher distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/F-distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FisherSnedecor : IDistribution
    {
        #region Private data
        // Distribution parameters
        private int d1;
        private int d2;

        // Derived values
        private double b;
        private double? mean;
        private double? variance;
        #endregion

        #region Fisher-Snedecor components
        /// <summary>
        /// Initializes the Fisher distribution.
        /// </summary>
        /// <param name="d1">First degree of freedom</param>
        /// <param name="d2">Second degree of freedom</param>
        public FisherSnedecor(int d1 = 1, int d2 = 1)
        {
            if (d1 <= 0)
                throw new ArgumentOutOfRangeException("d1", "The value must be greater than zero.");
            if (d2 <= 0)
                throw new ArgumentOutOfRangeException("d2", "The value must be greater than zero.");

            this.d1 = d1;
            this.d2 = d2;

            this.b = Special.Beta(d1 * 0.5, d2 * 0.5);
        }
        /// <summary>
        /// Gets the value of the first degree of freedom.
        /// </summary>
        public int D1
        {
            get { return d1; }
        }
        /// <summary>
        /// Gets the value of the second degree of freedom.
        /// </summary>
        public int D2
        {
            get { return d2; }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                if (!mean.HasValue)
                {
                    if (d2 <= 2)
                    {
                        mean = Double.NaN;
                    }
                    else
                    {
                        mean = d2 / (d2 - 2.0);
                    }
                }

                return mean.Value;
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
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                if (!variance.HasValue)
                {
                    if (d2 <= 4)
                    {
                        variance = Double.NaN;
                    }
                    else
                    {
                        variance = (2.0 * d2 * d2 * (d1 + d2 - 2)) /
                            (d1 * (d2 - 2) * (d2 - 2) * (d2 - 4));
                    }
                }

                return variance.Value;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (d1 > 2)
                {
                    double a = (d1 - 2.0) / d1;
                    double b = d2 / (d2 + 2.0);
                    return a * b;
                }

                return Double.NaN;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                if (d2 > 6)
                {
                    double v1 = 2 * d1 + d2 - 2;
                    double v2 = Math.Sqrt(8 * (d2 - 4));
                    double v3 = Math.Sqrt(d1 * (d1 + d2 - 2));
                    double v4 = d2 - 6;
                    return v1 * v2 / (v3 * v4);
                }
                return double.NaN;
            }
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
            get { return new RangeDouble(0, double.PositiveInfinity); }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x <= 0)
                return 0;

            double u = (d1 * x) / (d1 * x + d2);
            return Special.BetaIncomplete(d1 * 0.5, d2 * 0.5, u);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x <= 0)
                return 0;

            double u = Math.Pow(d1 * x, d1) * Math.Pow(d2, d2) /
                Math.Pow(d1 * x + d2, d1 + d2);
            return Math.Sqrt(u) / (x * b);
        }
        #endregion
    }
}
