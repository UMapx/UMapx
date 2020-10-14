using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the distribution of Nakagami.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Nakagami_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Nakagami : IDistribution
    {
        #region Private data
        private double mu;
        private double omega;
        private double constant; // 2 * μ ^ μ / (Γ(μ) * ω ^ μ))
        private double nratio;   // -μ / ω
        private double twoMu1;   // 2 * μ - 1.0
        #endregion

        #region Nakagami components
        /// <summary>
        ///Initializes the distribution of Nakagami.
        /// </summary>
        public Nakagami()
        {
            Initialize(0.5, 1);
        }
        /// <summary>
        /// Initializes the distribution of Nakagami.
        /// </summary>
        /// <param name="mu">Shape factor</param>
        /// <param name="omega">Spread rate</param>
        public Nakagami(double mu, double omega)
        {
            Initialize(mu, omega);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="mu"></param>
        /// <param name="omega"></param>
        private void Initialize(double mu, double omega)
        {
            Mu = mu;
            Omega = omega;

            double twoMuMu = 2.0 * Math.Pow(mu, mu);
            double gammaMu = Special.Gamma(mu);
            double spreadMu = Math.Pow(omega, mu);
            nratio = -mu / omega;
            twoMu1 = 2.0 * mu - 1.0;

            constant = twoMuMu / (gammaMu * spreadMu);
        }
        /// <summary>
        /// Gets or sets the value of the shape factor.
        /// </summary>
        public double Mu
        {
            get
            {
                return mu;
            }
            set
            {
                if (value <= 0.5)
                    throw new Exception("Invalid argument value");

                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the spread coefficient value.
        /// </summary>
        public double Omega
        {
            get
            {
                return omega;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.omega = value;
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
                return (Special.Gamma(mu + 0.5) / Special.Gamma(mu)) * Math.Sqrt(omega / mu);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                double a = Math.Sqrt(2) / 2;
                double b = ((2 * mu - 1) * omega) / mu;
                return a * Math.Sqrt(b);
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
                double a = Special.Gamma(mu + 0.5) / Special.Gamma(mu);
                return omega * (1.0 - (1.0 / mu) * (a * a));
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        ///  Gets the value of entropy.
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
            {
                return 0;
            }

            return Special.GammaP(mu, (mu / omega) * (x * x));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x <= 0)
            {
                return 0;
            }

            return constant * Math.Pow(x, twoMu1) * Math.Exp(nratio * x * x);
        }
        #endregion
    }
}
