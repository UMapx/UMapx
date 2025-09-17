using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the inverse chi-square distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Inverse-chi-squared_distribution
    /// </remarks>
    [Serializable]
    public class InverseChiSquare : IDistribution
    {
        #region Private data
        private int v = 1;
        #endregion

        #region Inverse chi-square components
        /// <summary>
        /// Initializes the inverse chi-square distribution.
        /// </summary>
        public InverseChiSquare() { }
        /// <summary>
        /// Initializes the inverse chi-square distribution.
        /// </summary>
        /// <param name="degreesOfFreedom">Degrees of freedom (positive integer)</param>
        public InverseChiSquare(int degreesOfFreedom)
        {
            V = degreesOfFreedom;
        }
        /// <summary>
        /// Gets or sets the degrees of freedom (positive integer).
        /// </summary>
        public int V
        {
            get { return this.v; }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.v = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get { return new RangeFloat(0f, float.PositiveInfinity); }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                if (v > 2)
                    return 1f / (v - 2f);

                return float.PositiveInfinity;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public float Variance
        {
            get
            {
                if (v <= 2)
                    return float.NaN;
                if (v <= 4)
                    return float.PositiveInfinity;

                float vf = v;
                float denom = (vf - 2f) * (vf - 2f) * (vf - 4f);
                return 2f / denom;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get { return float.NaN; }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get { return new float[] { 1f / (v + 2f) }; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public float Skewness
        {
            get { return float.NaN; }
        }
        /// <summary>
        /// Gets the excess kurtosis (kurtosis minus 3).
        /// </summary>
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get { return float.NaN; }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= 0f)
                return 0f;

            float vf = v;
            float a = Maths.Pow(2f, -vf / 2f);
            float b = Maths.Pow(x, -vf / 2f - 1f);
            float c = Maths.Exp(-1f / (2f * x));
            float d = Special.Gamma(vf / 2f);
            return (a * b * c) / d;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= 0f)
                return 0f;

            return Special.GammaQ(v / 2f, 1f / (2f * x));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Value</returns>
        public float Entropy
        {
            get
            {
                float vf = v;
                return vf / 2f + Maths.Log(0.5f) + Special.LogGamma(vf / 2f) - (1f - vf / 2f) * Special.DiGamma(vf / 2f);
            }
        }
        #endregion
    }
}
