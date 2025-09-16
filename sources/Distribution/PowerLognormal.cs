using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the power lognormal distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.itl.nist.gov/div898/handbook/eda/section3/eda366e.htm
    /// </remarks>
    [Serializable]
    public class PowerLognormal : IDistribution
    {
        #region Private data
        private float power;
        private float sigma;
        #endregion

        #region Constructor
        /// <summary>
        /// Initializes the power lognormal distribution with the given parameters.
        /// </summary>
        /// <param name="power">Shape parameter (must be greater than zero)</param>
        /// <param name="sigma">Scale parameter (must be greater than zero)</param>
        public PowerLognormal(float power, float sigma)
        {
            Power = power;
            Sigma = sigma;
        }
        #endregion

        #region Parameters
        /// <summary>
        /// Gets or sets the shape parameter (must be greater than zero).
        /// </summary>
        public float Power
        {
            get => power;
            set
            {
                if (value <= 0f)
                {
                    throw new ArgumentOutOfRangeException(nameof(Power), "Power must be greater than zero");
                }

                power = value;
            }
        }
        /// <summary>
        /// Gets or sets the scale parameter (must be greater than zero).
        /// </summary>
        public float Sigma
        {
            get => sigma;
            set
            {
                if (value <= 0f)
                {
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Sigma must be greater than zero");
                }

                sigma = value;
            }
        }
        #endregion

        #region Distribution properties
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support => new RangeFloat(0f, float.PositiveInfinity);
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        public float Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the value of differential entropy.
        /// </summary>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion

        #region Methods
        /// <summary>
        /// Returns the value of the probability density function f(x) = power · φ(z) · Φ(z)^(power - 1) / (x · sigma).
        /// Here z = ln(x) / sigma.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= 0f)
            {
                return 0f;
            }

            double z = Math.Log(x) / sigma;
            double pdf = StandardNormalPdf(z);
            double cdf = StandardNormalCdf(z);

            if (cdf <= double.Epsilon)
            {
                return 0f;
            }

            cdf = Math.Min(cdf, 1.0);
            double factor = power / (x * sigma);
            double value = factor * pdf * Math.Pow(cdf, power - 1.0);
            return (float)value;
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function F(x) = Φ(z)^power, where z = ln(x) / sigma.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= 0f)
            {
                return 0f;
            }

            double z = Math.Log(x) / sigma;
            double cdf = StandardNormalCdf(z);

            if (cdf <= 0.0)
            {
                return 0f;
            }

            if (cdf >= 1.0)
            {
                return 1f;
            }

            double result = Math.Pow(cdf, power);
            return (float)result;
        }
        #endregion

        #region Private helpers
        private static double StandardNormalPdf(double x)
        {
            const double invSqrt2Pi = 0.39894228040143267794;
            return invSqrt2Pi * Math.Exp(-0.5 * x * x);
        }

        private static double StandardNormalCdf(double x)
        {
            return 0.5 * (1.0 + Special.Erf((float)(x / Math.Sqrt(2.0))));
        }
        #endregion
    }
}
