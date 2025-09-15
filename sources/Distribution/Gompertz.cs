using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Gompertz distribution.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gompertz_distribution
    /// </remarks>
    [Serializable]
    public class Gompertz : IDistribution
    {
        #region Private data
        private float eta;
        private float b;
        #endregion

        #region Gompertz distribution
        /// <summary>
        ///Initializes the Gompertz distribution.
        /// </summary>
        /// <param name="eta">Form parameter η > 0</param>
        /// <param name="b">Scale parameter b > 0</param>
        public Gompertz(float eta, float b)
        {
            Eta = eta; B = b;
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter η > 0.
        /// </summary>
        public float Eta
        {
            get
            {
                return this.eta;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.eta = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter b > 0.
        /// </summary>
        public float B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.b = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                float e1 = -Special.Ei(-eta);
                return Maths.Exp(eta) * e1 / b;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        /// <remarks>
        /// Variance equals (2·e^{η}·E₁(η) − e^{η}·E₁(2η) − e^{η}·E₁(η)²) / b².
        /// </remarks>
        public float Variance
        {
            get
            {
                float e1 = -Special.Ei(-eta);
                float e2 = -Special.Ei(-2f * eta);
                float expEta = Maths.Exp(eta);
                return (2f * expEta * e1 - expEta * e2 - expEta * e1 * e1) / (b * b);
            }
        }
        /// <summary>
        /// Gets the mode values.
        /// </summary>
        public float[] Mode
        {
            get
            {
                if (eta >= 1)
                    return new float[] { 0f };

                return new float[] { (1 / b) * Maths.Log(1 / eta) };
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public float Median
        {
            get
            {
                return (1.0f / b) * Maths.Log((-1 / eta) * Maths.Log(0.5f) + 1);
            }
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
        /// <remarks>
        /// Full kurtosis equals 3 plus this value.
        /// </remarks>
        public float Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        public float Entropy
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
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x < 0)
            {
                return 0;
            }

            float ebx = Maths.Exp(b * x);
            return 1.0f - Maths.Exp(-eta * (ebx - 1.0f));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <remarks>
        /// The normalization ensures that the integral of the PDF over its support equals one.
        /// </remarks>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x < 0)
            {
                return 0;
            }

            float expBx = Maths.Exp(b * x);
            float exponent = b * x - eta * (expBx - 1f);
            return b * eta * Maths.Exp(exponent);
        }
        #endregion
    }
}
