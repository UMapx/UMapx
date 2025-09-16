using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the symmetric geometric (discrete Laplace) distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Laplace_distribution
    /// </remarks>
    [Serializable]
    public class SymmetricGeometric : IDistribution
    {
        #region Private data
        private float p = 0.5f;
        private float q = 0.5f;
        private float halfP = 0.25f;
        #endregion

        #region Symmetric geometric components
        /// <summary>
        /// Initializes the symmetric geometric distribution with p = 0.5.
        /// </summary>
        public SymmetricGeometric()
        {
        }
        /// <summary>
        /// Initializes the symmetric geometric distribution.
        /// </summary>
        /// <param name="p">Probability of zero (0, 1].</param>
        public SymmetricGeometric(float p)
        {
            P = p;
        }
        /// <summary>
        /// Gets or sets the probability assigned to zero (0, 1].
        /// </summary>
        public float P
        {
            get
            {
                return p;
            }
            set
            {
                if (value <= 0f || value > 1f)
                    throw new ArgumentException("Invalid argument value");

                p = value;
                q = 1f - p;
                halfP = 0.5f * p;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(float.NegativeInfinity, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value (always zero).
        /// </summary>
        public float Mean
        {
            get
            {
                return 0f;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (p == 0f)
                    return float.NaN;

                return ((2f - p) * (1f - p)) / (p * p);
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                return new float[] { 0f };
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return 0f;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get
            {
                return 0f;
            }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        public float Excess
        {
            get
            {
                if (q <= 0f)
                    return float.NaN;

                float numerator = -2f * p * p - 3f * p + 6f;
                float denominator = p * p - 3f * p + 2f;
                return numerator / denominator;
            }
        }
        /// <summary>
        /// Returns the value of the probability mass function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Function(float x)
        {
            int k = (int)Maths.Floor(x);
            if (x != k)
            {
                return 0f;
            }

            if (k == 0)
            {
                return p;
            }

            return halfP * Maths.Pow(q, Maths.Abs(k));
        }
        /// <summary>
        /// Returns the value of the cumulative distribution function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public float Distribution(float x)
        {
            int k = (int)Maths.Floor(x);

            if (k < 0)
            {
                int n = -k;
                return 0.5f * Maths.Pow(q, n);
            }

            return 1f - 0.5f * Maths.Pow(q, k + 1);
        }
        /// <summary>
        /// Gets the entropy value.
        /// </summary>
        public float Entropy
        {
            get
            {
                if (q == 0f)
                {
                    return 0f;
                }

                return -Maths.Log(p) + q * Maths.Log(2f) - (q / p) * Maths.Log(q);
            }
        }
        #endregion
    }
}
