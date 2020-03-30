// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2020
// Valery Asiryan
// Moscow, Russia

using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                             UMAPX.DISTRIBUTION
    // **************************************************************************
    // Designed by Valery Asiryan (c), 2015-2020
    // Moscow, Russia.
    // **************************************************************************

    #region Distributions
    /// <summary>
    /// Defines the Gaussian distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Normal_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Gaussian : IDistribution
    {
        #region Private data
        private double sigma = 1;
        private double mu = 0;
        #endregion;

        #region Gaussian components
        /// <summary>
        /// Initializes the Gaussian distribution.
        /// </summary>
        public Gaussian() { }
        /// <summary>
        /// Initializes the Gaussian distribution.
        /// </summary>
        /// <param name="sigma">Standard deviation</param>
        /// <param name="mu">Mathematical expectation</param>
        public Gaussian(double sigma, double mu)
        {
            Sigma = sigma;
            Mu = mu;
        }
        /// <summary>
        /// Gets or sets the standard deviation.
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Gets or sets the mathematical expectation.
        /// </summary>
        public double Mu
        {
            get
            {
                return this.mu;
            }
            set
            {
                this.mu = value;
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
            get
            {
                return this.mu;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return sigma * sigma;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return this.mu;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return this.mu;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0.0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 0.0;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return Maths.Exp(Maths.Pow((x - mu), 2) / (-2.0 * sigma * sigma)) / (Maths.Sqrt(2.0 * Maths.Pi) * sigma);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            return 0.5 + 0.5 * Special.Erf((x - mu) / Maths.Sqrt(2.0 * sigma * sigma));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(sigma * Maths.Sqrt(2 * Maths.Pi * Maths.E), Maths.E);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the logarithmic Gaussian distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Log-normal_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LogGaussian : IDistribution
    {
        #region Private data
        private double sigma = 1;
        private double mu = 0;
        #endregion;

        #region GaussianLog components
        /// <summary>
        /// Initializes the logarithmic Gaussian distribution.
        /// </summary>
        public LogGaussian() { }
        /// <summary>
        /// Initializes the logarithmic Gaussian distribution.
        /// </summary>
        /// <param name="sigma">Standard deviation</param>
        /// <param name="mu">Mathematical expectation</param>
        public LogGaussian(double sigma, double mu)
        {
            Sigma = sigma;
            Mu = mu;
        }
        /// <summary>
        /// Gets or sets the standard deviation.
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the mathematical expectation.
        /// </summary>
        public double Mu
        {
            get
            {
                return this.mu;
            }
            set
            {
                this.mu = value;
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
                return Maths.Exp(mu + sigma * sigma / 2);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return (Maths.Exp(sigma * sigma) - 1) * Maths.Exp(2 * mu + sigma * sigma);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Maths.Exp(mu);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return Maths.Exp(mu - sigma * sigma);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (Maths.Exp(sigma * sigma) + 2.0) * Maths.Sqrt(Maths.Exp(sigma * sigma) - 1.0);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return Maths.Exp(4 * sigma * sigma) + 2.0 * Maths.Exp(3 * sigma * sigma) + 3.0 * Maths.Exp(3 * sigma * sigma) - 6.0;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Maths.Exp(Maths.Pow((Maths.Log(x) - mu), 2) / (-2.0 * sigma * sigma)) / (Maths.Sqrt(2.0 * Maths.Pi) * sigma * x);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 0.5 + 0.5 * Special.Erf((Maths.Log(x) - mu) / Maths.Sqrt(sigma * 1.414));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return 0.5 + 0.5 * Maths.Log(2 * Maths.Pi * sigma * sigma) + mu;
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the Wiener semicircular distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Wigner : IDistribution
    {
        #region Private data
        private double r;
        #endregion;

        #region Wigner components
        /// <summary>
        /// Initializes the Wiener semicircular distribution.
        /// </summary>
        /// <param name="r">Radius</param>
        public Wigner(double r)
        {
            R = r;
        }
        /// <summary>
        /// Gets or sets the radius value.
        /// </summary>
        public double R
        {
            get
            {
                return this.r;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.r = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(-r, r);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
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
                return r * r / 4.0;
            }
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
                return 0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return -1;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (Math.Abs(x) > r)
            {
                return double.NaN;
            }

            double r2 = r * r, x2 = x * x;
            double a = Math.Sqrt(r2 - x2);
            double b = 2 / (Math.PI * r2);
            return b * a;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (Math.Abs(x) > r)
            {
                return double.NaN;
            }

            double r2 = r * r, x2 = x * x;
            double a = Math.Sqrt(r2 - x2);
            double b = x / (Math.PI * r2);
            double c = Math.Asin(x / r) / Maths.Pi;
            return 0.5 + b * a + c;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(Maths.Pi * r) - 1.0 / 2;
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the logarithmic distribution of Rayleigh.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Rayleigh_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Rayleigh : IDistribution
    {
        #region Private data
        private double sigma = 1;
        #endregion;

        #region Rayleigh components
        /// <summary>
        /// Initializes the Rayleigh logarithmic distribution.
        /// </summary>
        public Rayleigh() { }
        /// <summary>
        /// Initializes the Rayleigh logarithmic distribution.
        /// </summary>
        /// <param name="sigma">Scale parameter</param>
        public Rayleigh(double sigma)
        {
            Sigma = sigma;
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter.
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.sigma = value;
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
                return Maths.Sqrt(Maths.Pi / 2.0) * sigma;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return (2.0 - Maths.Pi / 2.0) * Maths.Pow(sigma);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return sigma * Maths.Sqrt(Maths.Log(4));
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return sigma;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2 * Maths.Sqrt(Maths.Pi) * (Maths.Pi - 3) / Maths.Pow(4 - Maths.Pi, 1.5);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return -(6 * Maths.Pi - 24 * Maths.Pi + 16) / Maths.Pow(4 - Maths.Pi);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return x / sigma / sigma * Maths.Exp(-(x * x) / (2 * sigma * sigma));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 1 - Maths.Exp(-(x * x) / (2 * sigma * sigma));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return 1 + Maths.Log(sigma / Maths.Log(2)) + Maths.Gamma / 2;
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the exponential distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Exponential_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Exponential : IDistribution
    {
        #region Private data
        private double l = 1;
        #endregion

        #region Exp components
        /// <summary>
        /// Initializes an exponential distribution.
        /// </summary>
        public Exponential() { }
        /// <summary>
        /// Initializes an exponential distribution.
        /// </summary>
        /// <param name="lambda">Intensity parameter (0, + inf)</param>
        public Exponential(double lambda)
        {
            Lambda = lambda;
        }
        /// <summary>
        /// Gets or sets the value of the intensity parameter (0, + inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.l = value;
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
                return Maths.Pow(l, -1);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return Maths.Pow(l, -2);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return 0.0;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Maths.Log(2) / l;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2.0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6.0;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return l * Maths.Exp(-l * x);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 1 - Maths.Exp(-l * x);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return 1 - Maths.Log(l);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the Cauchy distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Cauchy_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Cauchy : IDistribution
    {
        #region Private data
        private double g = 0.5;
        private double x0 = 0;
        #endregion

        #region Caushi components
        /// <summary>
        /// Initializes the Cauchy distribution.
        /// </summary>
        public Cauchy() { }
        /// <summary>
        /// Initializes the Cauchy distribution.
        /// </summary>
        /// <param name="gamma">Scale factor (0, + inf)</param>
        /// <param name="x0">Shift coefficient</param>
        public Cauchy(double gamma, double x0)
        {
            Gamma = gamma;
            X0 = x0;
        }
        /// <summary>
        /// Gets or sets the value of the scale factor (0, + inf).
        /// </summary>
        public double Gamma
        {
            get
            {
                return this.g;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.g = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the shift coefficient.
        /// </summary>
        public double X0
        {
            get
            {
                return this.x0;
            }
            set
            {
                this.x0 = value;
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
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return double.PositiveInfinity;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return x0;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return x0;
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
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return 1.0 / (Maths.Pi * g * (1.0 + Maths.Pow((x - x0) / g)));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            return 1.0 / Maths.Pi * Maths.Atg((x - x0) / g) + 0.5;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(4 * Maths.Pi * g);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the Weibull distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Weibull_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Weibull : IDistribution
    {
        #region Private data
        private double l = 1;
        private double k = 1;
        #endregion

        #region Weibull components
        /// <summary>
        /// Initializes the Weibull distribution.
        /// </summary>
        public Weibull() { }
        /// <summary>
        /// Initializes the Weibull distribution.
        /// </summary>
        /// <param name="lambda">Scale factor (0, + inf)</param>
        /// <param name="k">Shape factor (0, + inf)</param>
        public Weibull(double lambda, double k)
        {
            Lambda = lambda;
        }
        /// <summary>
        /// Gets or sets the value of the scale factor (0, + inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.l = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the form factor (0, + inf).
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                this.k = Maths.Max(0.0000001, value);
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
                return l * Special.Gamma(1.0 + 1.0 / k);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return l * l * (Special.Gamma(1.0 + 2.0 / k) - Special.Gamma(1.0 + 1.0 / k));
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return l * Maths.Pow(k - 1, 1.0 / k) / Maths.Pow(k, 1.0 / k);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return l * Maths.Pow(Maths.Log(2), 1.0 / k);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (Special.Gamma(1.0 + 3.0 / k) * Maths.Pow(l, 3) - 3 * Mean * Special.Gamma(1 + 2.0 / k) * Maths.Pow(l, 2) + 2 * Maths.Pow(Mean, 3)) / (Variance * Maths.Sqrt(Variance));
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return (Special.Gamma(1.0 + 4.0 / k) * Maths.Pow(l, 4) - 4 * Mean * Special.Gamma(1 + 3.0 / k) * Maths.Pow(l, 3) + 6 * Maths.Pow(Mean, 2) * Math.Pow(l, 2) * Special.Gamma(1.0 + 2.0 / k) - 3 * Maths.Pow(Mean, 4)) / (Variance * Variance);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return (k / l) * Maths.Pow(x / l, k - 1) * Maths.Exp(-Maths.Pow(x / l, k));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 1 - Maths.Exp(-Maths.Pow(x / l, k));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Maths.Gamma * (1.0 - 1.0 / k) + Maths.Pow(l / k, k) + Maths.Log(l / k);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the Laplace distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Laplace_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Laplace : IDistribution
    {
        #region Private data
        private double a = 1;
        private double b = 0;
        #endregion

        #region Laplace components
        /// <summary>
        /// Initializes the Laplace distribution.
        /// </summary>
        public Laplace() { }
        /// <summary>
        /// Initializes the Laplace distribution.
        /// </summary>
        /// <param name="alfa">Scale factor (0, + inf)</param>
        /// <param name="beta">Shift coefficient</param>
        public Laplace(double alfa, double beta)
        {
            Alfa = alfa;
            Beta = beta;
        }
        /// <summary>
        /// Gets or sets the value of the scale factor (0, + inf).
        /// </summary>
        public double Alfa
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the shift coefficient.
        /// </summary>
        public double Beta
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
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
            get
            {
                return b;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return 2.0 / a / a;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return b;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return b;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 3;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return a / 2.0 * Maths.Exp(-a * Maths.Abs(x - b));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x <= b)
            {
                return 0.5 * Maths.Exp(a * (x - b));
            }
            return 1 - 0.5 * Maths.Exp(-a * (x - b));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(4 * Maths.Pi * a);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the geometric distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Geometric_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Geometric : IDistribution
    {
        #region Private data
        private double p = 0.2;
        private double q = 0.8;
        #endregion

        #region Geometric components
        /// <summary>
        /// Initializes the geometric distribution.
        /// </summary>
        public Geometric() { }
        /// <summary>
        /// Initializes the geometric distribution.
        /// </summary>
        /// <param name="p">Probability of "success" [0, 1]</param>
        public Geometric(double p)
        {
            P = p;
        }
        /// <summary>
        /// Gets or sets the probability value of "success" [0, 1].
        /// </summary>
        public double P
        {
            get
            {
                return this.p;
            }
            set
            {
                if (value < 0 || value > 1)
                    throw new ArgumentException("Invalid argument value");

                this.p = value;
                this.q = 1 - this.p;
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
                return q / p;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return q / p / p;
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
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (2 - p) / Maths.Sqrt(1 - p);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6 + p * p / q;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Maths.Pow(q, x) * p;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 1 - Maths.Pow(q, x);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return -Maths.Log2(p) - q / p * Maths.Log2(q);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the Poisson distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Poisson_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Poisson : IDistribution
    {
        #region Private data
        private double l = 1;
        #endregion

        #region Poisson components
        /// <summary>
        /// Initializes the Poisson distribution.
        /// </summary>
        public Poisson() { }
        /// <summary>
        /// Initializes the Poisson distribution.
        /// </summary>
        /// <param name="lambda">Parameter λ (0, +inf)</param>
        public Poisson(double lambda)
        {
            Lambda = lambda;
        }
        /// <summary>
        /// Gets or sets the value of the parameter λ (0, +inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.l = value;
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
                return l;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return l;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return Math.Floor(l);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Floor(l + 0.333 - 0.02 / l);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return Math.Pow(l, -0.5);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return Math.Pow(l, -1.0);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Math.Exp(-l) * Math.Pow(l, x) / Special.Factorial(x);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Special.GammaQ(x + 1, l);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return l * (1 - Math.Log(l)) + Math.Exp(-l) * Row(l);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="l"></param>
        /// <returns></returns>
        private double Row(double l)
        {
            double sum = 0;
            int k, n = 20;
            double fac;

            for (k = 0; k < n; k++)
            {
                fac = Special.Factorial(k);
                sum += Math.Pow(l, k) * Math.Log(fac) / fac;
            }

            return sum;
        }
        #endregion
    }
    /// <summary>
    /// Defines the uniform distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Uniform : IDistribution
    {
        #region Private data
        private double a = 1;
        private double b = 1;
        #endregion

        #region Uniform components
        /// <summary>
        /// Initializes the uniform distribution.
        /// </summary>
        public Uniform() { }
        /// <summary>
        /// Initializes the uniform distribution.
        /// </summary>
        /// <param name="a">Shift parameter a</param>
        /// <param name="b">Shift parameter b</param>
        public Uniform(double a, double b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the shift parameter a.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the shift parameter b.
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
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
            get
            {
                return (a + b) / 2.0;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return Math.Pow(b - a, 2)  / 12.0;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return (a - b) / 2.0;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Mean;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return -1.2;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < a)
            {
                return 0;
            }
            else if (x > b)
            {
                return 0;
            }
            return 1.0 / (b - a);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < a)
            {
                return 0;
            }
            else if (x > b)
            {
                return 1;
            }
            return (x - a) / (b - a);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return Math.Log(b - a);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the beta distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Beta_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Beta : IDistribution
    {
        #region Private data
        private double a = 1;
        private double b = 1;
        #endregion

        #region Beta components
        /// <summary>
        /// Initializes beta distribution.
        /// </summary>
        public Beta() { }
        /// <summary>
        /// Initializes beta distribution.
        /// </summary>
        /// <param name="a">Parameter a</param>
        /// <param name="b">Parameter b</param>
        public Beta(double a, double b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the parameter a.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the parameter b.
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, 1);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return a / (a + b);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return (a * b) / Maths.Pow(a + b) / (a + b + 1);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (a > 1 && b > 1)
                {
                    return (a - 1) / (a + b - 2);
                }
                return double.NaN;
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
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2 * (b - a) * Math.Sqrt(a + b + 1) / (a + b + 2) / Math.Sqrt(a * b);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                double a2 = a * a, b2 = b * b, a3 = a2 * a;
                return 6 * (a3 - a2 * (2 * b - 1) + b2 * (b + 1) - 2 * a * b * (b + 2)) / (a * b * (a + b + 2) * (a + b + 3));
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x > 1)
            {
                return 0;
            }
            else if (x < 0)
            {
                return 0;
            }
            return Math.Pow(x, a - 1) * Math.Pow(1 - x, b - 1) / Special.Beta(a, b);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x > 1)
            {
                return 0;
            }
            else if (x < 0)
            {
                return 0;
            }
            return Special.Beta(a, b, x);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
    /// <summary>
    /// Defines the Gamma-distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gamma_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Gamma : IDistribution
    {
        #region Private data
        private double thetta = 1;
        private double k = 1;
        #endregion

        #region Gamma components
        /// <summary>
        /// Initializes the Gamma-distribution.
        /// </summary>
        public Gamma() { }
        /// <summary>
        /// Initializes the Gamma-distribution.
        /// </summary>
        /// <param name="thetta">Parameter θ (0, +inf)</param>
        /// <param name="k">Parameter k (0, +inf)</param>
        public Gamma(double thetta, double k)
        {
            Thetta = thetta; K = k;
        }
        /// <summary>
        /// Gets or sets the parameter θ (0, +inf).
        /// </summary>
        public double Thetta
        {
            get
            {
                return this.thetta;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.thetta = value;
            }
        }
        /// <summary>
        /// Gets or sets the parameter k (0, +inf).
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.k = value;
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
                return thetta * k;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return k * thetta * thetta;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (k >= 1)
                {
                    return (k - 1) * thetta;
                }
                return double.NaN;
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
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2 / Math.Sqrt(k);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6.0 / k;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Math.Pow(x, k - 1) * Math.Exp(-x / thetta) / (Special.Gamma(k) * Math.Pow(thetta, k));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Special.GammaP(k, x / thetta);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return k * thetta + (1 - k) * Math.Log(thetta) + Special.GammaLog(k); // + (1 - k) * Special.Ksi(k);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the Pareto distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Pareto_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Pareto : IDistribution
    {
        #region Private data
        private double xm = 1;
        private double k = 1;
        #endregion

        #region Pareto components
        /// <summary>
        /// Initializes the Pareto distribution.
        /// </summary>
        public Pareto() { }
        /// <summary>
        /// Initializes the Pareto distribution.
        /// </summary>
        /// <param name="xm">Scale factor θ (0, +inf)</param>
        /// <param name="k">Parameter k (0, +inf)</param>
        public Pareto(double xm, double k)
        {
            Xm = xm; K = k;
        }
        /// <summary>
        /// Gets or sets the scale factor Xm (0, +inf).
        /// </summary>
        public double Xm
        {
            get
            {
                return this.xm;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.xm = value;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor k (0, +inf).
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.k = value;
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
                return (k * xm) / (k - 1);
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                if (k > 2)
                {
                    return Maths.Pow(xm / k - 1) * (k / (k - 2));
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return xm;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return xm * Maths.Sqrt(2, k);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                if (k > 3)
                {
                    return 2 * (1 + k) / (k - 3) * Math.Sqrt((k - 2) / k);
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
                double k2 = k * k;
                double k3 = k2 * k;
                return 6 * (k3 + k2 + 6 * k - 2) / (k * (k - 3) * (k - 4));
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < xm)
            {
                return 0;
            }
            return k * Math.Pow(xm, k) / Math.Pow(x, k + 1);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < xm)
            {
                return 0;
            }
            return 1.0 - Math.Pow(xm / x, k);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return k * xm + (1 - k) * Math.Log(xm) + Special.GammaLog(k); // + (1 - k) * Special.Ksi(k);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the distribution of Bernoulli.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Bernoulli_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Bernoulli : IDistribution
    {
        #region Private data
        private double p = 0.2;
        private double q = 0.8;
        #endregion

        #region Bernoully components
        /// <summary>
        /// Initializes a Bernoulli distribution.
        /// </summary>
        public Bernoulli() { }
        /// <summary>
        /// Initializes a Bernoulli distribution.
        /// </summary>
        /// <param name="p">Probability of success [0, 1]</param>
        public Bernoulli(double p)
        {
            P = p;
        }
        /// <summary>
        /// Gets or sets the probability of success [0, 1].
        /// </summary>
        public double P
        {
            get
            {
                return this.p;
            }
            set
            {
                if (value < 0 || value > 1)
                    throw new ArgumentException("Invalid argument value");

                this.p = value;
                this.q = 1 - this.p;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, 1);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                return p;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return p * q;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (q > p)
                {
                    return 0;
                }
                else if (q < p)
                {
                    return 1;
                }
                return 0.5;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                if (q > p)
                {
                    return 0;
                }
                else if (q < p)
                {
                    return 1;
                }
                return 0.5;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (p - q) / Maths.Sqrt(p * q);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return (6 * p * p - 6 * p + 1) / p * (1 - p);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x == 0)
            {
                return q;
            }
            return p;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            else if (x >= 1)
            {
                return 1;
            }
            return q;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return -q * Maths.Log(q) - p * Maths.Log(p);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the log-logistic distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Log-logistic_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LogLogistic : IDistribution
    {
        #region Private data
        private double a = 1;
        private double b = 1;
        #endregion

        #region LogLogistic components
        /// <summary>
        /// Initializes the log-logistic distribution.
        /// </summary>
        public LogLogistic() { }
        /// <summary>
        /// Initializes the log-logistic distribution.
        /// </summary>
        /// <param name="a">Parameter a</param>
        /// <param name="b">Parameter b</param>
        public LogLogistic(double a, double b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the value of parameter a.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of parameter b.
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.b = value;
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
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (b > 1)
                {
                    return a * Math.Pow((b - 1) / (b + 1), 1 / b);
                }
                return 0;
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
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return (b / a) * Math.Pow(x / a, b - 1) / (1.0 + Math.Pow(Math.Pow(x / a, b), 2));
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            return 1.0 / (1 + Math.Pow(x / a, -b));
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
    /// <summary>
    /// Defines the binomial distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Binomial_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Binomial : IDistribution
    {
        #region Private data
        private double n = 20;
        private double p = 0.5;
        private double q = 0.5;
        #endregion

        #region Binomial components
        /// <summary>
        /// Initializes the binomial distribution.
        /// </summary>
        public Binomial() { }
        /// <summary>
        /// Initializes the binomial distribution.
        /// </summary>
        /// <param name="n">Number of experiments (>0)</param>
        /// <param name="p">Probability of success [0, 1]</param>
        public Binomial(double n, double p)
        {
            N = n; P = p;
        }
        /// <summary>
        /// Gets or sets number of experiments.
        /// </summary>
        public double N
        {
            get 
            { 
                return n; 
            }
            set
            {
                this.n = value;
            }
        }
        /// <summary>
        /// Gets or sets probability of success [0, 1].
        /// </summary>
        public double P
        {
            get 
            { 
                return p; 
            }
            set
            {
                if (value > 1 || value < 0)
                    throw new Exception("Invalid argument value");

                this.p = value;
                this.q = 1.0 - p;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, this.n);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get 
            { 
                return n * p; 
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return n * p * q;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                double test = (n + 1) * p;

                if (test <= 0 || (int)test != test)
                    return Math.Floor(test);

                if (test <= n)
                    return test;

                if (test == n + 1)
                    return n;

                return Double.NaN;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Floor(n * p);
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (q - p) / Math.Sqrt(n * p * q);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return (1 - 6 * p * q) / Variance;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0 || x > n)
            {
                return 0;
            }

            double a = Special.LogBinomial(n, x);
            double b = x == 0 ? 0 : x * Math.Log(p);
            double c = (n - x);
            double d = Math.Log(1 - p);
            double log = a + b + c * d;

            return Math.Exp(log);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
                return 0;
            if (x >= n)
                return 1;

            double a = n - x;
            double b = x + 1;
            return Special.Beta(a, b, q);
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
    /// <summary>
    /// Defines the hypergeometric distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hypergeometric_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Hypergeometric : IDistribution
    {
        #region Private data
        private double n = 30;
        private double k = 20;
        private double d = 20;
        #endregion

        #region Hypergeometric components
        /// <summary>
        /// Initializes the hypergeometric distribution.
        /// </summary>
        public Hypergeometric() { }
        /// <summary>
        /// Initializes the hypergeometric distribution.
        /// </summary>
        /// <param name="n">Parameter N [0, +inf]</param>
        /// <param name="k">Parameter D [0, N]</param>
        /// <param name="d">Parameter K [0, N]</param>
        public Hypergeometric(double n, double k, double d) 
        {
            N = n; K = k; D = d;
        }
        /// <summary>
        /// Gets or sets the value of the parameter N [0, +inf].
        /// </summary>
        public double N
        {
            get
            {
                return n;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.n = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter D [0, N].
        /// </summary>
        public double D
        {
            get
            {
                return d;
            }
            set
            {
                if (value < 0 || value > N)
                    throw new Exception("Invalid argument value");

                this.d = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter k [0, N].
        /// </summary>
        public double K
        {
            get
            {
                return k;
            }
            set
            {
                if (value < 0 || value > N)
                    throw new Exception("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, k);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get 
            { 
                return k * d / n; 
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return k * (d / n) * ((n - d) / n) * ((n - k) / (n - 1.0)); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                double num = (k + 1) * (d + 1);
                double den = n + 2;
                return Math.Floor(num / den);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                double k1 = (n - 2 * d) * Math.Pow(n - 1, 0.5) * (n - 2 * k);
                double k2 = (n * d * (n - d) * Math.Pow(n - k, 0.5) * (n - 2));
                return k1 / k2;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                double n2 = n * n;
                double k1 = (n2 * (n - 1)) / (k * (n - 2) * (n - 3) * (n - k));
                double k2 = (n * (n + 1) - 6 * n * (n - k)) / (d * (n - d));
                double k3 = (3 * n * (n - k) * (n + 6)) / n2 - 6;
                return k1 * (k2 + k3);
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < Math.Max(0, k + d - n) || x > Math.Min(d, k))
            {
                return 0;
            }

            double a = Special.Binomial(d, x);
            double b = Special.Binomial(n - d, k - x);
            double c = Special.Binomial(n, k);
            return (a * b) / c;
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            double sum = 0;
            int k = (int)x;
            for (int i = 0; i <= k; i++)
            {
                sum += Function(i);
            }

            return sum;
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
    /// <summary>
    /// Defines the logistic distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Logistic_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Logistic : IDistribution
    {
        #region Private data
        private double mu = 5;
        private double s = 2;
        #endregion

        #region Logistic components
        /// <summary>
        /// Initializes the logistic distribution.
        /// </summary>
        /// <param name="mu">Parameter μ</param>
        /// <param name="s">Parameter s (0, +inf]</param>
        public Logistic(double mu, double s)
        {
            Mu = mu; S = s;
        }
        /// <summary>
        /// Initializes the logistic distribution.
        /// </summary>
        public Logistic() { }
        /// <summary>
        /// Gets or sets the value of the parameter μ.
        /// </summary>
        public double Mu
        { 
            get 
            { 
                return mu; 
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter s (0, +inf].
        /// </summary>
        public double S
        {
            get 
            {
                return s; 
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.s = value;
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
            get 
            { 
                return mu; 
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return (s * s * Math.PI * Math.PI) / 3.0; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return mu;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return mu; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 1.2;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                return Math.Log(s) + 2;
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            double z = (x - mu) / s;
            return 1.0 / (1 + Math.Exp(-z));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double z = (x - mu) / s;
            double num = Math.Exp(-z);
            double a = (1 + num);
            double den = s * a * a;

            return num / den;
        }
        #endregion
    }
    /// <summary>
    /// Defines the Rademacher distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Rademacher_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Rademacher : IDistribution
    {
        #region Rademacher components
        /// <summary>
        /// Initializes the Rademacher distribution.
        /// </summary>
        public Rademacher() { }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(-1, 1);
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
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return 1; }
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
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return double.NaN; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return -2;
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < -1)
            {
                return 0;
            }
            if (x >= 1)
            {
                return 1;
            }
            return 0.5;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x == -1)
            {
                return 0.5;
            }
            if (x == 1)
            {
                return 0.5;
            }
            return 0.0;
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                return Math.Log(2);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the triangular distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Triangular_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Triangular : IDistribution
    {
        #region Private data
        private double a;
        private double b;
        private double c;
        #endregion

        #region Triangular components
        /// <summary>
        /// Initializes the triangular distribution.
        /// </summary>
        public Triangular() { }
        /// <summary>
        /// Initializes the triangular distribution.
        /// </summary>
        /// <param name="a">Parameter a ∈ (-inf, +inf)</param>
        /// <param name="b">Parameter b ∈ (-inf, +inf)</param>
        /// <param name="c">Parameter c ∈ (-inf, +inf)</param>
        public Triangular(double a, double b, double c)
        {
            A = a; B = b; C = c;
        }
        /// <summary>
        /// Gets or sets the value of the parameter a ∈ (-inf, +inf).
        /// </summary>
        public double A
        { 
            get 
            { 
                return a; 
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter b ∈ (-inf, +inf).
        /// </summary>
        public double B 
        { 
            get 
            { 
                return b; 
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter c ∈ (-inf, +inf).
        /// </summary>
        public double C
        {
            get
            {
                return c;
            }
            set
            {
                this.c = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(a, b);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return (a + b + c) / 3.0; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return (a * a + b * b + c * c - a * b - a * c - b * c) / 18; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                double median;
                if (c >= (a + b) / 2.0)
                {
                    median = a + Math.Sqrt((b - a) * (c - a)) / 1.4142135623730950488016887242097;
                }
                else
                {
                    median = b - Math.Sqrt((b - a) * (b - c)) / 1.4142135623730950488016887242097;
                }

                return median;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return c; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                double k1 = (a + b - 2 * c) * (2 * a - b - c) * (a - 2 * b + c);
                double k2 = 5 * (a * a + b * b + c * c - a * b - a * c - b * c);
                return 1.4142135623730950488016887242097 * k1 / Math.Pow(k2, 1.5);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return -0.6;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                return 0.5 + Math.Log((b - a) / 2);
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < a)
                return 0;

            if (x >= a && x <= c)
                return ((x - a) * (x - a)) / ((b - a) * (c - a));

            if (x > c && x <= b)
                return 1 - ((b - x) * (b - x)) / ((b - a) * (b - c));

            return 1;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < a)
                return 0;

            if (x >= a && x <= c)
                return (2 * (x - a)) / ((b - a) * (c - a));

            if (x > c && x <= b)
                return (2 * (b - x)) / ((b - a) * (b - c));

            return 0;
        }
        #endregion
    }
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
    /// <summary>
    /// Defines the Levy distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/L%C3%A9vy_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Levy : IDistribution
    {
        #region Prviate data
        private double mu = 0;
        private double scale = 1;
        #endregion

        #region Levy components
        /// <summary>
        /// Initializes the Levy distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ</param>
        /// <param name="c">Scale factor (>0)</param>
        public Levy(double mu, double c)
        {
            Mu = mu; C = c;
        }
        /// <summary>
        /// Gets or sets the shift factor.
        /// </summary>
        public double Mu
        {
            get 
            { 
                return mu; 
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor (> 0).
        /// </summary>
        public double C
        {
            get 
            { 
                return scale; 
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.scale = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(mu, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return Double.PositiveInfinity; }
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
            get { return Double.PositiveInfinity; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (mu == 0)
                {
                    return scale / 3.0;
                }
                return Double.NaN;
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
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                return (1.0 + 3.0 * Maths.Gamma + Math.Log(16 * Math.PI * scale * scale)) / 2.0;
            }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < mu)
            {
                return 0;
            }

            return Special.Erfc(Math.Sqrt(scale / (2 * (x - mu))));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < mu)
            {
                return 0;
            }
            double z = x - mu;
            double a = Math.Sqrt(scale / (2.0 * Math.PI));
            double b = Math.Exp(-(scale / (2 * z)));
            double c = Math.Pow(z, 3.0 / 2.0);

            return a * b / c;
        }
        #endregion
    }
    /// <summary>
    /// Defines the logarithmic distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Logarithmic_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Logarithmic : IDistribution
    {
        #region Private data
        private double p = 0.66;
        #endregion

        #region Logarithmic components
        /// <summary>
        /// Initializes the logarithmic distribution.
        /// </summary>
        /// <param name="p">Parameter</param>
        public Logarithmic(double p)
        {
            P = p;
        }
        /// <summary>
        /// Gets or sets the value of the parameter p ∈ (0, 1].
        /// </summary>
        public double P
        {
            get 
            { 
                return p; 
            }
            set 
            {
                if (value <= 0 || value > 1)
                    throw new Exception("Invalid argument value");

                this.p = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(1, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get 
            { 
                return (-1 / Math.Log(1 - p)) * p / (1 - p); 
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get 
            {
                double k1 = p + Math.Log(1 - p);
                double k2 = Math.Pow(1 - p, 2) * Math.Pow(Math.Log(1 - p), 2);
                return -p * k1 / k2;
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
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return 1; }
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
            if (x > 1)
            {
                return 0;
            }
            return 1 + Special.Beta(x + 1, 0) / Math.Log(1 - p);
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
            if (x > 1)
            {
                return 0;
            }
            return -1 / Math.Log(1 - p) * Math.Pow(p, x) / x;
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion
    }
    /// <summary>
    /// Defines the beta distribution of the second kind.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Beta_prime_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class BetaPrime : IDistribution
    {
        #region Private data
        private double alpha = 1; // shape (α)
        private double beta = 1;  // shape (β)
        #endregion

        #region Beta-prime components
        /// <summary>
        /// Initializes beta distribution of the second kind.
        /// </summary>
        /// <param name="alpha">Parameter α (0, +inf)</param>
        /// <param name="beta">Parameter β (0, +inf)</param>
        public BetaPrime(double alpha, double beta)
        {
            Alpha = alpha; Beta = beta;
        }
        /// <summary>
        /// Gets or sets the value of the parameter α ∈ (0, +inf).
        /// </summary>
        public double Alpha
        {
            get 
            { 
                return alpha; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.alpha = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter β ∈ (0, +inf).
        /// </summary>
        public double Beta
        {
            get 
            { 
                return beta; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.beta = value;
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
                if (beta > 1)
                {
                    return alpha / (beta - 1);
                }

                return Double.PositiveInfinity;
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
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (alpha >= 1)
                {
                    return (alpha - 1) / (beta + 1);
                }

                return 0.0;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                if (beta > 2.0)
                {
                    double num = alpha * (alpha + beta - 1);
                    double den = (beta - 2) * Math.Pow(beta - 1, 2);
                    return num / den;
                }
                else if (beta > 1.0)
                {
                    return Double.PositiveInfinity;
                }

                return double.NaN;
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
            {
                return 0;
            }
            return Special.BetaIncomplete(alpha, beta, x / (1 + x));
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

            double num = Math.Pow(x, alpha - 1) * Math.Pow(1 + x, -alpha - beta);
            double den = Special.Beta(alpha, beta);
            return num / den;
        }
        #endregion
    }
    /// <summary>
    /// Defines the Birnbaum-Saunders distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Birnbaum–Saunders_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class BirnbaumSaunders : IDistribution
    {
        #region Private data
        private double mu = 0;
        private double beta = 1;
        private double gamma = 1;
        #endregion

        #region Birnbaum-Saunders components
        /// <summary>
        /// Initializes the Birnbaum-Saunders distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ ∈ (0, +inf)</param>
        /// <param name="beta">Scale factor β ∈ (0, +inf).</param>
        /// <param name="gamma">Shape factor γ ∈ (0, +inf)</param>
        public BirnbaumSaunders(double mu, double beta, double gamma)
        {
            Mu = mu; Beta = beta; Gamma = gamma;
        }
        /// <summary>
        /// Gets or sets the shift factor μ ∈ (0, +inf).
        /// </summary>
        public double Mu
        {
            get 
            { 
                return mu; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor β ∈ (0, +inf).
        /// </summary>
        public double Beta
        {
            get 
            { 
                return beta; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.beta = value;
            }
        }
        /// <summary>
        /// Gets or sets the form factor γ ∈ (0, +inf).
        /// </summary>
        public double Gamma
        {
            get 
            { 
                return gamma; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.gamma = value;
            }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support 
        {
            get
            {
                return new RangeDouble(this.mu, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get 
            { 
                return 1 + 0.5 * gamma * gamma; 
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
            get { return gamma * gamma * (1 + (5 * gamma * gamma) / 4); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
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
            double a = Math.Sqrt(x);
            double b = Math.Sqrt(1.0 / x);
            double z = (a - b) / gamma;

            // Normal cumulative distribution function
            return Special.Erfc(-z / 1.4142135623731) * 0.5;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double c = x - mu;

            double a = Math.Sqrt(c / beta);
            double b = Math.Sqrt(beta / c);

            double alpha = (a + b) / (2 * gamma * c);
            double z = (a - b) / gamma;

            // Normal cumulative distribution function
            return alpha * Special.Erfc(-z / 1.4142135623731) * 0.5;
        }
        #endregion
    }
    /// <summary>
    /// Defines the xi-square distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Chi-squared_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ChiSquare : IDistribution
    {
        #region Private data
        private int k = 1;
        #endregion

        #region Chi-square components
        /// <summary>
        /// Initializes the xi-square distribution.
        /// </summary>
        /// <param name="k">Degrees of freedom (0, +inf)</param>
        public ChiSquare(int k)
        {
            K = k;
        }
        /// <summary>
        /// Gets or sets the degrees of freedom (0, +inf).
        /// </summary>
        public int K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.k = value;
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
                return k;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return k - 2.0 / 3;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (k >= 2)
                {
                    return k - 2;
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return 2 * k;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return Math.Sqrt(8.0 / k);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 12.0 / k;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get
            {
                double k2 = k / 2.0;
                double s1 = Math.Log(2.0 * Special.Gamma(k2));
                double s2 = (1.0 - k2) * Special.DiGamma(k2);
                return k2 + s1 + s2;
            }
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
            return Special.GammaP(k / 2.0, x / 2.0);
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
                return 1;
            }
            return Special.GammaQ(k / 2.0, x / 2.0);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Gumbel distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gumbel_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Gumbel : IDistribution
    {
        #region Private data
        private double mu = 0;
        private double beta = 1;
        #endregion

        #region Gumbel components
        /// <summary>
        /// Initializes the Gumbel distribution.
        /// </summary>
        /// <param name="mu">Shear rate μ ∈ (-inf, +inf)</param>
        /// <param name="beta">Scale factor β ∈ (0, +inf).</param>
        public Gumbel(double mu, double beta)
        {
            Mu = mu; Beta = beta;
        }
        /// <summary>
        /// Gets or sets the shift factor μ ∈ (-inf, +inf).
        /// </summary>
        public double Mu
        {
            get
            {
                return mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor β ∈ (0, +inf).
        /// </summary>
        public double Beta
        {
            get
            {
                return beta;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.beta = value;
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
            get
            {
                return mu + beta * Maths.Gamma;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return mu - beta * Math.Log(Math.Log(2));
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return ((Math.PI * Math.PI) / 6.0) * beta * beta; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return mu;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 1.14;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 12.0 / 5;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { return Math.Log(beta) + Maths.Gamma + 1; }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            double z = (x - mu) / beta;
            return Math.Exp(-Math.Exp(-z));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double z = (x - mu) / beta;
            return (1 / beta) * Math.Exp(-(z + Math.Exp(-z)));
        }
        #endregion
    }
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
    /// <summary>
    /// Defines the U-quadratic distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/U-quadratic_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class UQuadratic : IDistribution
    {
        #region Private data
        double a;
        double b;
        double alpha;
        double beta;
        #endregion

        #region UQuadratic components
        /// <summary>
        /// Initializes the U-quadratic distribution.
        /// </summary>
        /// <param name="a">Parameter a ∈ (0, +inf)</param>
        /// <param name="b">Parameter b ∈ (a, +inf)</param>
        public UQuadratic(double a, double b)
        {
            A = a; B = b;
            this.alpha = 12 / Math.Pow(b - a, 3);
            this.beta = (b + a) / 2;
        }
        /// <summary>
        /// Gets or sets the value of the parameter a ∈ (0, +inf).
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter b ∈ (a, +inf).
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (b < a)
                    throw new Exception("The value of parameter b must be either greater than or equal to a");

                this.b = value;
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return (int)a & (int)b;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return (a + b) / 2.0d; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { return (a + b) / 2.0d; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return (3.0d / 20.0) * Math.Pow(b - a, 2.0); }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 3.0 / 112.0 * Math.Pow(b - a, 4.0);
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(a, b); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < a)
                return 0;

            if (x > b)
                return 1;

            return (alpha / 3) * (Math.Pow(x - beta, 3) + Math.Pow(beta - a, 3));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < a)
                return 0;

            if (x > b)
                return 0;

            return alpha * Math.Pow(x - beta, 2);
        }
        #endregion
    }
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
    /// <summary>
    /// Defines the distribution of Erlang.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Erlang_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Erlang : IDistribution
    {
        #region Private data
        private int k = 1;
        private double lambda = 0.5;
        #endregion

        #region Erlang distribution
        /// <summary>
        /// Initializes the distribution of Erlang.
        /// </summary>
        /// <param name="k">Form parameter k ∈ (0, +inf)</param>
        /// <param name="lambda">λ-parameter λ ∈ (0, +inf)</param>
        public Erlang(int k, double lambda)
        {
            K = k; Lambda = lambda;
        }
        /// <summary>
        /// Gets or sets the value of the parameter k ∈ (0, +inf).
        /// </summary>
        public int K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter λ ∈ (0, +inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.lambda;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.lambda = value;
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
                return this.k / this.lambda;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return double.NaN;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                return this.k / Math.Pow(this.lambda, 2.0);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return (this.k - 1) / this.lambda;
            }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2.0 / Math.Sqrt(k);
            }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6.0 / k;
            }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        public double Entropy
        {
            get
            {
                return (1 - k) * Special.DiGamma(k) + Math.Log(Special.Gamma(k) / lambda) + k;
            }
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return double.NaN;
            }
            return Math.Pow(lambda, k) * Math.Pow(x, k - 1) * Math.Exp(-lambda * x) / Special.Factorial(k - 1);
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return double.NaN;
            }
            return Special.GammaIncomplete(k, lambda * x) / Special.Factorial(k - 1);
        }
        #endregion
    }
    /// <summary>
    /// Defines the compact Cauchy distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wrapped_Cauchy_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class WrappedCauchy : IDistribution
    {
        #region Private data
        private double mu;
        private double gamma;
        #endregion

        #region Wrapped distribution
        /// <summary>
        /// Initializes the compact Cauchy distribution.
        /// </summary>
        /// <param name="mu">Parameter μ</param>
        /// <param name="gamma">Parameter γ > 0</param>
        public WrappedCauchy(double mu, double gamma)
        {
            this.mu = mu;
            this.gamma = gamma;
        }
        /// <summary>
        /// Gets or sets the value of the parameter μ.
        /// </summary>
        public double Mu
        {
            get
            {
                return this.mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the parameter γ > 0.
        /// </summary>
        public double Gamma
        {
            get
            {
                return this.gamma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.gamma = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return mu; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return 1 - Math.Exp(-gamma); }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(-Math.PI, Math.PI); }
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
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { return Math.Log(2 * Math.PI * (1 - Math.Exp(-2 * gamma))); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double constant = (1.0 / (2 * Math.PI));
            return constant * Math.Sinh(gamma) / (Math.Cosh(gamma) - Math.Cos(x - mu));
        }
        #endregion
    }
    /// <summary>
    /// Defines the distribution of Kumaraswa.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Kumaraswamy_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Kumaraswamy : IDistribution
    {
        #region Private data
        private double a;
        private double b;
        #endregion

        #region Distribution
        /// <summary>
        /// Initializes the distribution of Kumarasva.
        /// </summary>
        /// <param name="a">Form parameter a > 0</param>
        /// <param name="b">Form parameter b > 0</param>
        public Kumaraswamy(double a, double b)
        {
            this.A = a;
            this.B = b;
        }
        /// <summary>
        /// Gets or sets form parameter a > 0.
        /// </summary>
        public double A
        {
            get
            {
                return a;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.a = value;
            }
        }
        /// <summary>
        /// Gets or sets form parameter b > 0.
        /// </summary>
        public double B
        {
            get
            {
                return b;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.b = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get
            {
                double num = b * Special.Gamma(1 + (1 / a)) * Special.Gamma(b);
                double den = Special.Gamma(1 + (1 / a) + b);

                return num / den;
            }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get
            {
                double alpha = momentGeneratingFunction(2, a, b);
                double beta = Math.Pow(momentGeneratingFunction(1, a, b), 2);
                return alpha - beta;
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Pow(1 - Math.Pow(2, -1 / b), 1 / a);
            }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {

                if ((a >= 1) && (b >= 1) && (a != 1 && b != 1))
                {
                    double num = a - 1;
                    double den = a * b - 1;
                    return Math.Pow(num / den, 1 / a);
                }

                return Double.NaN;
            }
        }
        /// <summary>
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { return double.NaN; }
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
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(0, 1); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x > 1)
                return 1;

            if (x < 0)
                return 0;

            double xa = Math.Pow(x, a);
            return 1 - Math.Pow(1 - xa, b);
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x > 1)
                return 0;

            if (x < 0)
                return 0;

            return a * b * Math.Pow(x, a - 1) * Math.Pow(1 - Math.Pow(x, a), b - 1);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        protected static double momentGeneratingFunction(int n, double a, double b)
        {
            return (b * Special.Beta(1.0 + ((double)n) / a, b));
        }
        #endregion
    }
    /// <summary>
    /// Defines the Gompertz distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gompertz_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Gompertz : IDistribution
    {
        #region Private data
        private double eta;
        private double b;
        #endregion

        #region Gompertz distribution
        /// <summary>
        ///Initializes the Gompertz distribution.
        /// </summary>
        /// <param name="eta">Form parameter η > 0</param>
        /// <param name="b">Scale parameter b > 0</param>
        public Gompertz(double eta, double b)
        {
            Eta = eta; B = b;
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter η > 0.
        /// </summary>
        public double Eta
        {
            get
            {
                return this.eta;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.eta = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter b > 0.
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.b = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                if (eta >= 1)
                    return 0;

                return (1 / b) * Math.Log(1 / eta);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return (1.0 / b) * Math.Log((-1 / eta) * Math.Log(0.5) + 1);
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
        /// Returns the value of differential entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(0, Double.PositiveInfinity); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }

            double ebx = Math.Exp(b * x);
            return 1.0 - Math.Exp(-eta * (ebx - 1.0));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }

            double a1 = b * eta * Math.Exp(eta);
            double a2 = Math.Exp(b * x);
            double a3 = Math.Exp(-eta * Math.Exp(b * x));
            return a1 * a2 * a3;
        }
        #endregion
    }
    /// <summary>
    /// Defines the hyperbolic secant distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class HyperbolicSecant : IDistribution
    {
        #region Distribution
        /// <summary>
        /// Initializes the hyperbolic secant distribution.
        /// </summary>
        public HyperbolicSecant() { }
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
            get { return 0; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return 1; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get { return 2; }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(Double.NegativeInfinity, Double.PositiveInfinity); }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        public double Entropy
        {
            get { return (4.0 / Math.PI) * Maths.G; }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            double angle = Math.Atan(Math.Exp(x * Math.PI / 2.0));
            return 2 * angle / Math.PI;
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return 0.5 * Maths.Sch(x * (Math.PI / 2.0));
        }
        #endregion
    }
    /// <summary>
    /// Defines the arcsine distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Arcsine_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Arcsine : IDistribution
    {
        #region Distribution
        /// <summary>
        /// Initializes the arcsine distribution.
        /// </summary>
        public Arcsine() { }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return 0.5; }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { return 0.5; }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { return 1.0 / 8; }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { return double.NaN; }
        }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        public double Skewness
        {
            get { return 0; }
        }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        public double Excess
        {
            get { return -1.5; }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(0, 1); }
        }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        public double Entropy
        {
            get { return Math.Log(Math.PI / 4); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            return 2.0 / Math.PI * Math.Asin(Math.Sqrt(x));
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            if (x < 1 && x > 0)
            {
                return 1.0 / (Math.PI * Math.Sqrt(x * (1 - x)));
            }
            return double.NaN;
        }
        #endregion
    }
    /// <summary>
    /// Defines the Burr distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Burr_distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Burr : IDistribution
    {
        #region Private data
        private double c;
        private double k;
        #endregion

        #region Burr distribution
        /// <summary>
        /// Initializes the Burr distribution.
        /// </summary>
        /// <param name="c">Form parameter c > 0</param>
        /// <param name="k">Scale parameter k > 0</param>
        public Burr(double c, double k)
        {
            C = c; K = k;
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter c > 0.
        /// </summary>
        public double C
        {
            get
            {
                return this.c;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.c = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the scale parameter k > 0.
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.k = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { return k * Special.Beta(k - 1.0 / c, 1.0 + 1.0 / c); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get
            {
                return Math.Pow((c - 1) / (k * c + 1), 1.0 / c);
            }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Pow(Math.Pow(2, 1.0 / k) - 1.0, 1.0 / c);
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
        /// Returns the value of differential entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(double.Epsilon, Double.PositiveInfinity); }
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
                return double.NaN;
            }

            return 1.0 - Math.Pow(1.0 + Math.Pow(x, c), -k);
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
                return double.NaN;
            }

            double a = c * k;
            double b = Math.Pow(x, c - 1);
            double d = 1 + Math.Pow(x, c);
            return a * b / Math.Pow(d, k + 1);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Fisher's Z-distribution.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Fisher%27s_z-distribution
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FisherZ : IDistribution
    {
        #region Private data
        private double d1;
        private double d2;
        #endregion

        #region FisherZ distribution
        /// <summary>
        /// Initializes the Fisher Z-distribution.
        /// </summary>
        /// <param name="d1">Degree of freedom d1 > 0</param>
        /// <param name="d2">Degree of freedom d2 > 0</param>
        public FisherZ(double d1, double d2)
        {
            D1 = d1; D2 = d2;
        }
        /// <summary>
        /// Gets or sets the degree of freedom d1 > 0.
        /// </summary>
        public double D1
        {
            get
            {
                return this.d1;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.d1 = value;
            }
        }
        /// <summary>
        /// Gets or sets the degree of freedom d2 > 0.
        /// </summary>
        public double D2
        {
            get
            {
                return this.d2;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.d2 = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
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
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
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
        /// Returns the value of differential entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(Double.NegativeInfinity, Double.PositiveInfinity); }
        }
        /// <summary>
        /// Returns the value of the probability distribution function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the probability density function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            // helpers:

            double d12 = d1 / 2.0;
            double d22 = d2 / 2.0;

            // first equation:

            double a = Math.Pow(d1, d12);
            double b = Math.Pow(d2, d22);
            double c = 2 * a * b;
            double d = Special.Beta(d12, d22);
            double e = c / d;

            // second equation:

            double f = Math.Exp(d1 * x);
            double g = d1 * Math.Exp(2 * x) + d2;
            double h = d12 + d22;
            double j = f / Math.Pow(g, h);

            // result of F(x, d1, d2):
            return e * j;
        }
        #endregion
    }
    #endregion

    #region Unsorted distributions
    /// <summary>
    /// Defines the distribution of the conical shape.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Cone-shape_distribution_function
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ConeShape : IDistribution
    {
        #region Private data
        private double a;
        #endregion

        #region Cone-Shape components
        /// <summary>
        /// Initializes the distribution of the conical shape.
        /// </summary>
        /// <param name="a">Coefficient</param>
        public ConeShape(double a = 0.001)
        {
            A = a;
        }
        /// <summary>
        /// Gets or sets the coefficient value.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
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
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the kernel density function.
        /// </summary>
        /// <param name="eta">Argument</param>
        /// <param name="tau">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double eta, double tau)
        {
            double ksi = Maths.Pi * eta * tau;
            double psi = Math.Exp(-2 * Maths.Pi * a * tau * tau);
            return Math.Sin(ksi) / ksi * psi;
        }
        /// <summary>
        /// Returns the value of the kernel distribution function.
        /// </summary>
        /// <param name="t">Argument</param>
        /// <param name="tau">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double t, double tau)
        {
            if (Math.Abs(tau) >= 2 * Math.Abs(t))
            {
                return 1.0 / tau * Math.Exp(-2 * Maths.Pi * a * tau * tau);
            }
            return 0;
        }
        #endregion
    }
    /// <summary>
    /// Defines the distribution of Choi Williams.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Choi%E2%80%93Williams_distribution_function
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ChoiWilliams : IDistribution
    {
        #region Private data
        private double a;
        #endregion

        #region Choi-Williams components
        /// <summary>
        /// Initializes the Choi-Williams distribution.
        /// </summary>
        /// <param name="a">Coefficient</param>
        public ChoiWilliams(double a = 0.001)
        {
            A = a;
        }
        /// <summary>
        /// Gets or sets the coefficient value.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
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
        /// Gets the value of entropy.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Returns the value of the kernel density function.
        /// </summary>
        /// <param name="eta">Argument</param>
        /// <param name="tau">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double eta, double tau)
        {
            double ksi = eta * tau;
            return Math.Exp(-a * ksi * ksi);
        }
        /// <summary>
        /// Returns the value of the kernel distribution function.
        /// </summary>
        /// <param name="t">Argument</param>
        /// <param name="tau">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Distribution(double t, double tau)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    #endregion

    #region Probabilities
    /// <summary>
    /// Defines the Bayes probability class.
    /// </summary>
    [Serializable]
    public class Bayes
    {
        #region Private data
        private double[] Pp;
        private double Pa;
        private int N;
        #endregion

        #region Bayes components
        /// <summary>
        /// Initializes the Bayes probability class.
        /// </summary>
        /// <param name="stat">Array of statistical probabilities</param>
        /// <param name="prior">An array of a priori probabilities (before experiment)</param>
        public Bayes(double[] stat, double[] prior)
        {
            if (stat.Length != prior.Length)
                throw new Exception("Arrays must be of the same dimensions.");


            this.N = prior.Length;
            this.Pp = new double[N];
            this.Pa = 0;
            int i;

            for (i = 0; i < N; i++)
            {
                Pa += stat[i] * prior[i];
            }

            for (i = 0; i < N; i++)
            {
                Pp[i] = stat[i] * prior[i] / Pa;
            }
        }
        /// <summary>
        /// Returns the value of the total probability.
        /// </summary>
        public double General
        {
            get
            {
                return Pa;
            }
        }
        /// <summary>
        /// Returns an array of values of posterior probabilities (after the experiment).
        /// </summary>
        public double[] Probabilities
        {
            get
            {
                return Pp;
            }
        }
        #endregion
    }
    #endregion

    #region Interface
    /// <summary>
    /// Defines the distribution interface.
    /// </summary>
    public interface IDistribution
    {
        #region Components
        /// <summary>
        /// Gets the support interval of the argument.
        /// </summary>
        RangeDouble Support { get; }
        /// <summary>
        /// Gets the mean value.
        /// </summary>
        double Mean { get; }
        /// <summary>
        /// Gets the variance value.
        /// </summary>
        double Variance { get; }
        /// <summary>
        /// Gets the median value.
        /// </summary>
        double Median { get; }
        /// <summary>
        /// Gets the mode value.
        /// </summary>
        double Mode { get; }
        /// <summary>
        /// Gets the value of the asymmetry coefficient.
        /// </summary>
        double Skewness { get; }
        /// <summary>
        /// Gets the kurtosis coefficient.
        /// </summary>
        double Excess { get; }
        /// <summary>
        /// Returns the value of differential entropy.
        /// </summary>
        /// <returns>Double precision floating point number</returns>
        double Entropy { get; }
        #endregion
    }
    #endregion
}
