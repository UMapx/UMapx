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
        private float b;
        private float? mean;
        private float? variance;
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

            this.b = Special.Beta(d1 * 0.5f, d2 * 0.5f);
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
        public float Mean
        {
            get
            {
                if (!mean.HasValue)
                {
                    if (d2 <= 2)
                    {
                        mean = float.NaN;
                    }
                    else
                    {
                        mean = d2 / (d2 - 2.0f);
                    }
                }

                return mean.Value;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (!variance.HasValue)
                {
                    if (d2 <= 4)
                    {
                        variance = float.NaN;
                    }
                    else
                    {
                        variance = (2.0f * d2 * d2 * (d1 + d2 - 2)) /
                            (d1 * (d2 - 2) * (d2 - 2) * (d2 - 4));
                    }
                }

                return variance.Value;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                if (d1 > 2)
                {
                    float a = (d1 - 2.0f) / d1;
                    float b = d2 / (d2 + 2.0f);
                    return a * b;
                }

                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                if (d2 > 6)
                {
                    float v1 = 2 * d1 + d2 - 2;
                    float v2 = (float)Math.Sqrt(8 * (d2 - 4));
                    float v3 = (float)Math.Sqrt(d1 * (d1 + d2 - 2));
                    float v4 = d2 - 6;
                    return v1 * v2 / (v3 * v4);
                }
                return float.NaN;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(0, float.PositiveInfinity); }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>float precision floating point number</returns>
        public float Distribution(float x)
        {
            if (x <= 0)
                return 0;

            float u = (d1 * x) / (d1 * x + d2);
            return Special.BetaIncomplete(d1 * 0.5f, d2 * 0.5f, u);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>float precision floating point number</returns>
        public float Function(float x)
        {
            if (x <= 0)
                return 0;

            float u = (float)Math.Pow(d1 * x, d1) * (float)Math.Pow(d2, d2) /
                (float)Math.Pow(d1 * x + d2, d1 + d2);
            return (float)Math.Sqrt(u) / (x * b);
        }
        #endregion
    }
}
