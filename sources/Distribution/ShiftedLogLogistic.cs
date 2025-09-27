using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the shifted log-logistic distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Shifted_log-logistic_distribution
    /// </remarks>
    [Serializable]
    public class ShiftedLogLogistic : IDistribution
    {
        #region Private data
        private float mu;
        private float sigma;
        private float ksi;
        #endregion

        #region Constructor
        /// <summary>
        /// Initializes the shifted log-logistic distribution with zero location, unit scale and zero shape.
        /// </summary>
        public ShiftedLogLogistic()
        {
            this.mu = 0f;
            this.sigma = 1f;
            this.ksi = 0f;
        }
        /// <summary>
        /// Initializes the shifted log-logistic distribution with the given location parameter.
        /// </summary>
        /// <param name="location">Location parameter</param>
        public ShiftedLogLogistic(float location)
        {
            this.mu = location;
            this.sigma = 1f;
            this.ksi = 0f;
        }
        /// <summary>
        /// Initializes the shifted log-logistic distribution with the given location and scale parameters.
        /// </summary>
        /// <param name="location">Location parameter</param>
        /// <param name="scale">Scale parameter (must be greater than zero)</param>
        public ShiftedLogLogistic(float location, float scale)
        {
            this.mu = location;
            this.sigma = scale;
            this.ksi = 0f;
        }
        /// <summary>
        /// Initializes the shifted log-logistic distribution with the given parameters.
        /// </summary>
        /// <param name="location">Location parameter</param>
        /// <param name="scale">Scale parameter (must be greater than zero)</param>
        /// <param name="shape">Shape parameter</param>
        public ShiftedLogLogistic(float location, float scale, float shape)
        {
            this.mu = location;
            this.sigma = scale;
            this.ksi = shape;
        }
        #endregion

        #region Parameters
        /// <summary>
        /// Gets or sets the location parameter.
        /// </summary>
        public float Mu
        {
            get => mu;
            set => mu = value;
        }
        /// <summary>
        /// Gets or sets the scale parameter (must be greater than zero).
        /// </summary>
        public float Sigma
        {
            get => sigma;
            set
            {
                if (value <= 0)
                {
                    throw new ArgumentOutOfRangeException(nameof(sigma), "Scale must be positive");
                }

                sigma = value;
            }
        }
        /// <summary>
        /// Gets or sets the shape parameter.
        /// </summary>
        public float Ksi
        {
            get => ksi;
            set => ksi = value;
        }
        #endregion

        #region Distribution properties
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                if (ksi > 0f)
                {
                    return new RangeFloat(mu - sigma / ksi, float.PositiveInfinity);
                }

                if (ksi < 0f)
                {
                    return new RangeFloat(float.NegativeInfinity, mu - sigma / ksi);
                }

                return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                if (Maths.Abs(ksi) < 1e-6f)
                {
                    return mu;
                }

                double s = sigma;
                double shape = ksi;
                double b = Math.PI * shape;
                double sin = Math.Sin(b);

                if (Math.Abs(sin) < 1e-12)
                {
                    return float.NaN;
                }

                double c = s / shape;
                double mean = mu + c * (b / sin - 1.0);
                return (float)mean;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (Maths.Abs(ksi) < 1e-6f)
                {
                    return (float)(sigma * sigma * Math.PI * Math.PI / 3.0);
                }

                double shape = ksi;
                double pb = Math.PI * shape;
                double sinPb = Math.Sin(pb);
                double sin2Pb = Math.Sin(2.0 * pb);

                if (Math.Abs(sinPb) < 1e-12 || Math.Abs(sin2Pb) < 1e-12)
                {
                    return float.NaN;
                }

                double a = sigma * sigma / (shape * shape);
                double b = 2.0 * pb / sin2Pb;
                double c = pb / sinPb;
                double variance = a * (b - c * c);
                return (float)variance;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median => mu;
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                if (Maths.Abs(ksi) < 1e-6f)
                {
                    return new float[] { mu };
                }

                double shape = ksi;
                if (Math.Abs(1.0 + shape) < 1e-8 || Math.Abs(1.0 - shape) < 1e-8)
                {
                    return new float[] { float.NaN };
                }

                double ratio = (1.0 - shape) / (1.0 + shape);
                if (ratio <= 0)
                {
                    return new float[] { float.NaN };
                }

                double pow = Math.Pow(ratio, shape);
                double mode = mu + (sigma / shape) * (pow - 1.0);
                return new float[] { (float)mode };
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness => float.NaN;
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        public float Excess => float.NaN;
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
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (Maths.Abs(ksi) < 1e-6f)
            {
                float z = (x - mu) / sigma;
                float e = Maths.Exp(-z);
                float denom = sigma * (1f + e);
                return e / (denom * (1f + e));
            }
            else
            {
                float zeta = (x - mu) / sigma;
                float baseValue = 1f + ksi * zeta;

                if ((ksi > 0f && baseValue <= 0f) || (ksi < 0f && baseValue <= 0f))
                {
                    return 0f;
                }

                double power = Math.Pow(baseValue, -1.0 / ksi - 1.0);
                double denom = sigma * Math.Pow(1.0 + Math.Pow(baseValue, -1.0 / ksi), 2.0);
                return (float)(power / denom);
            }
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (Maths.Abs(ksi) < 1e-6f)
            {
                float z = (x - mu) / sigma;
                return 1f / (1f + Maths.Exp(-z));
            }

            float zeta = (x - mu) / sigma;
            float baseValue = 1f + ksi * zeta;

            if (ksi > 0f && baseValue <= 0f)
            {
                return 0f;
            }

            if (ksi < 0f && baseValue <= 0f)
            {
                return 1f;
            }

            double value = 1.0 / (1.0 + Math.Pow(baseValue, -1.0 / ksi));
            return (float)value;
        }
        #endregion
    }
}
