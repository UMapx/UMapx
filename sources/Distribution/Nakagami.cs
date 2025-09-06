using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the distribution of Nakagami.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Nakagami_distribution
    /// </remarks>
    [Serializable]
    public class Nakagami : IDistribution
    {
        #region Private data
        private float mu;
        private float omega;
        private float constant; // 2 * μ ^ μ / (Γ(μ) * ω ^ μ))
        private float nratio;   // -μ / ω
        private float twoMu1;   // 2 * μ - 1.0
        #endregion

        #region Nakagami components
        /// <summary>
        ///Initializes the distribution of Nakagami.
        /// </summary>
        public Nakagami()
        {
            Initialize(0.5f, 1);
        }
        /// <summary>
        /// Initializes the distribution of Nakagami.
        /// </summary>
        /// <param name="mu">Shape factor</param>
        /// <param name="omega">Spread rate</param>
        public Nakagami(float mu, float omega)
        {
            Initialize(mu, omega);
        }
        /// <summary>
        /// Configures distribution parameters and precomputes constants.
        /// </summary>
        /// <param name="mu">Shape factor</param>
        /// <param name="omega">Spread coefficient</param>
        private void Initialize(float mu, float omega)
        {
            Mu = mu;
            Omega = omega;

            float twoMuMu = 2.0f * Maths.Pow(mu, mu);
            float gammaMu = Special.Gamma(mu);
            float spreadMu = Maths.Pow(omega, mu);
            nratio = -mu / omega;
            twoMu1 = 2.0f * mu - 1.0f;

            constant = twoMuMu / (gammaMu * spreadMu);
        }
        /// <summary>
        /// Gets or sets the value of the shape factor.
        /// </summary>
        public float Mu
        {
            get
            {
                return mu;
            }
            set
            {
                if (value <= 0.5)
                    throw new ArgumentException("Invalid argument value");

                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the spread coefficient value.
        /// </summary>
        public float Omega
        {
            get
            {
                return omega;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Invalid argument value");

                this.omega = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeFloat Support
        {
            get
            {
                return new RangeFloat(0, float.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public float Mean
        {
            get
            {
                return Special.Gamma(mu + 0.5f) / Special.Gamma(mu) * Maths.Sqrt(omega / mu);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public float Mode
        {
            get
            {
                float a = Maths.Sqrt2 / 2;
                float b = ((2 * mu - 1) * omega) / mu;
                return a * Maths.Sqrt(b);
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
                float a = Special.Gamma(mu + 0.5f) / Special.Gamma(mu);
                return omega * (1.0f - (1.0f / mu) * (a * a));
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
        /// Gets the kurtosis coefficient.
        /// </summary>
        public float Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        ///  Gets the value of entropy.
        /// </summary>
        public float Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Distribution(float x)
        {
            if (x <= 0)
            {
                return 0;
            }

            return Special.GammaP(mu, (mu / omega) * (x * x));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            if (x <= 0)
            {
                return 0;
            }

            return constant * Maths.Pow(x, twoMu1) * Maths.Exp(nratio * x * x);
        }
        #endregion
    }
}
