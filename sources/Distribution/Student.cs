using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Student's distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Student%27s_t-distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Student : IDistribution
    {
        #region Private data
        private double lambda;
        private double degrees;
        #endregion

        #region Student components
        /// <summary>
        /// Initializes the Student's distribution.
        /// </summary>
        /// <param name="n">Degrees of freedom n ∈ (0, +inf)</param>
        public Student(double n)
        {
            this.N = n;
            double num = Special.GammaLog((n + 1) / 2.0);
            double den = 0.5 * Math.Log(n * Maths.Pi) + Special.GammaLog(n / 2.0);
            this.lambda = num - den;
        }
        /// <summary>
        /// Gets or sets degrees of freedom n ∈ (0, +inf).
        /// </summary>
        public double N
        {
            get
            {
                return this.degrees;
            }
            set
            {
                if (value < 1)
                    throw new Exception("Invalid argument value");

                this.degrees = value;
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
            get { return 0; }
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
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                if (degrees > 2)
                    return degrees / (degrees - 2);
                else if (degrees > 1)
                    return Double.PositiveInfinity;
                return Double.NaN;
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
                if (degrees > 3)
                {
                    return 0;
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                if (degrees > 4)
                {
                    return 6.0 / (degrees - 4);
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                double a = Special.DiGamma((1 + degrees) / 2.0);
                double b = Special.DiGamma(degrees / 2.0);
                double c = (degrees + 1) / 2.0 * (a - b);
                double d = Math.Sqrt(degrees) * Special.Beta(degrees / 2.0, 0.5);

                return c + Math.Log(d);
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            double v = degrees;
            double sqrt = Math.Sqrt(x * x + v);
            double u = (x + sqrt) / (2 * sqrt);
            return Special.BetaIncomplete(v / 2.0, v / 2.0, u);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return Math.Exp(LogFunction(x));
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private double LogFunction(double x)
        {
            return lambda - ((degrees + 1) / 2.0) * Math.Log((x * x) / degrees);
        }
        #endregion
    }
}
