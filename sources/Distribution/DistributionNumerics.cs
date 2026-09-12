using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    /// <summary>
    /// Provides shared double-precision probability kernels and transformed quadrature.
    /// </summary>
    internal static class DistributionNumerics
    {
        /// <summary>Euler–Mascheroni constant used by entropy formulas</summary>
        internal const double EulerGamma = 0.57721566490153286061;
        /// <summary>Natural logarithm of sqrt(2*pi), used in Gaussian normalizations</summary>
        internal const double LogSqrtTwoPi = 0.91893853320467274178;

        /// <summary>
        /// Evaluates log(1 + x) with compensated addition near zero.
        /// </summary>
        /// <param name="x">Real argument greater than -1; -1 gives the logarithmic pole.</param>
        /// <returns>The natural logarithm of 1 + x.</returns>
        internal static double Log1p(double x)
        {
            double sum = 1 + x;
            if (sum == 1) return x;
            if (double.IsPositiveInfinity(sum)) return sum;
            return Math.Log(sum) * x / (sum - 1);
        }

        /// <summary>
        /// Evaluates log(1 + exp(x)) without overflowing the positive exponential.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>The nonnegative value log(1 + exp(x)).</returns>
        internal static double Softplus(double x) => x > 0 ? x + Log1p(Math.Exp(-x)) : Log1p(Math.Exp(x));

        /// <summary>
        /// Returns log(Phi(-x)), using a Mills-ratio expansion before the Gaussian tail underflows.
        /// </summary>
        /// <param name="x">Standard-normal threshold.</param>
        /// <returns>The logarithm of the standard-normal upper-tail probability.</returns>
        internal static double LogNormalSurvival(double x)
        {
            if (x < 0) return Log1p(-0.5 * Special.DistributionErfc(-x / Math.Sqrt(2)));
            if (x < 26) return Math.Log(0.5 * Special.DistributionErfc(x / Math.Sqrt(2)));
            return -0.5 * x * x + LogNormalMillsRatio(x);
        }

        /// <summary>
        /// Evaluates the logarithm of the scaled normal tail exp(x^2/2)*Phi(-x) by its decreasing asymptotic terms.
        /// </summary>
        /// <param name="x">Positive threshold, at least 26 at the call sites.</param>
        /// <returns>log(Phi(-x)) + x^2/2, retaining the Gaussian normalization constant.</returns>
        private static double LogNormalMillsRatio(double x)
        {
            // Mills ratio expansion, evaluated before the Gaussian tail underflows.
            double term = 1, sum = 1, square = x * x;
            for (int i = 1; i < 50; i++)
            {
                double next = -term * (2 * i - 1) / square;
                if (Math.Abs(next) >= Math.Abs(term)) break;
                sum += next; term = next;
                if (Math.Abs(term) < Math.Abs(sum) * 2e-16) break;
            }
            return -Math.Log(x) - LogSqrtTwoPi + Math.Log(sum);
        }

        /// <summary>
        /// Evaluates power*phi(x)*Phi(-x)^(power-1), combining Gaussian exponents to preserve extreme tails.
        /// </summary>
        /// <param name="x">Standardized real observation.</param>
        /// <param name="power">Positive power-normal shape parameter.</param>
        /// <returns>The power-normal probability density at x.</returns>
        internal static double PowerNormalDensity(double x, double power)
        {
            if (double.IsInfinity(x)) return 0;
            double logTail = LogNormalSurvival(x);
            // Combine the Gaussian exponents first; this avoids infinity minus infinity for small powers.
            double logDensity = x >= 26 ? -0.5 * power * x * x + Math.Log(power) - LogSqrtTwoPi
                + (power - 1) * LogNormalMillsRatio(x) :
                Math.Log(power) - 0.5 * x * x - LogSqrtTwoPi + (power - 1) * logTail;
            return Math.Exp(logDensity);
        }

        /// <summary>
        /// Returns beta-prime excess kurtosis, with positive infinity when the fourth moment diverges.
        /// </summary>
        /// <param name="a">Positive first shape parameter.</param>
        /// <param name="b">Positive second shape parameter.</param>
        /// <returns>Excess kurtosis, or positive infinity when b is at most four.</returns>
        internal static double BetaPrimeExcess(double a, double b)
        {
            if (b <= 4) return double.PositiveInfinity;
            return 6 * (a * (a + b - 1) * (5 * b - 11) + (b - 1) * (b - 1) * (b - 2)) /
                (a * (a + b - 1) * (b - 3) * (b - 4));
        }

        /// <summary>
        /// Returns the logarithm of the Poisson probability at a nonnegative integer count.
        /// </summary>
        /// <param name="lambda">Positive Poisson rate.</param>
        /// <param name="k">Nonnegative integer count.</param>
        /// <returns>log(P(X = k)) for a Poisson variable with mean lambda.</returns>
        internal static double PoissonLogMass(double lambda, double k)
        {
            return k == 0 ? -lambda : -lambda + k * Math.Log(lambda) - Special.DistributionLogGamma(k + 1);
        }

        // Adaptive integration on a finite transformed interval. The callers choose decaying tails.
        /// <summary>
        /// Integrates a real function on a finite interval with adaptive Simpson refinement.
        /// </summary>
        /// <param name="f">Integrand evaluated at real quadrature nodes.</param>
        /// <param name="a">Left integration endpoint.</param>
        /// <param name="b">Right integration endpoint.</param>
        /// <param name="tolerance">Absolute quadrature error budget.</param>
        /// <returns>The estimated integral over [a, b].</returns>
        internal static double Integrate(Func<double, double> f, double a, double b, double tolerance = 1e-11)
        {
            double mid = (a + b) / 2, fa = f(a), fm = f(mid), fb = f(b);
            return Refine(f, a, b, fa, fm, fb, (b - a) * (fa + 4 * fm + fb) / 6, tolerance, 24);
        }

        /// <summary>
        /// Refines a Simpson panel, reusing sampled values and splitting the absolute error budget between subintervals.
        /// </summary>
        /// <param name="f">Integrand evaluated at real quadrature nodes.</param>
        /// <param name="a">Left endpoint of the current panel.</param>
        /// <param name="b">Right endpoint of the current panel.</param>
        /// <param name="fa">Integrand value at the left endpoint.</param>
        /// <param name="fm">Integrand value at the midpoint.</param>
        /// <param name="fb">Integrand value at the right endpoint.</param>
        /// <param name="whole">Simpson estimate for the unsplit interval.</param>
        /// <param name="tolerance">Absolute quadrature error budget.</param>
        /// <param name="depth">Remaining subdivision depth.</param>
        /// <returns>The corrected integral estimate for the current panel.</returns>
        private static double Refine(Func<double, double> f, double a, double b, double fa, double fm, double fb,
            double whole, double tolerance, int depth)
        {
            double mid = (a + b) / 2, fl = f((a + mid) / 2), fr = f((mid + b) / 2);
            double left = (mid - a) * (fa + 4 * fl + fm) / 6;
            double right = (b - mid) * (fm + 4 * fr + fb) / 6;
            double error = left + right - whole;
            if (double.IsNaN(error)) return double.NaN;
            if (depth == 0 || Math.Abs(error) <= 15 * tolerance) return left + right + error / 15;
            return Refine(f, a, mid, fa, fl, fm, left, tolerance / 2, depth - 1)
                 + Refine(f, mid, b, fm, fr, fb, right, tolerance / 2, depth - 1);
        }

        /// <summary>
        /// Computes the Gompertz mean or variance through log(Exp(1)) quadrature, centering the variance integrand.
        /// </summary>
        /// <remarks>Uses X = log(1 + E/eta)/rate for E distributed exponentially with unit rate. The transformed variable log(E) is integrated over [-100, 4], split at zero to resolve the central mass.</remarks>
        /// <param name="eta">Positive Gompertz shape parameter.</param>
        /// <param name="rate">Positive Gompertz rate parameter.</param>
        /// <param name="variance">True selects variance; false selects the mean.</param>
        /// <returns>The requested moment in the original rate scale.</returns>
        internal static double GompertzMoment(double eta, double rate, bool variance)
        {
            double logEta = Math.Log(eta), scale = Math.Max(1, eta);
            Func<double, double> value = t => scale * Softplus(t - logEta);
            Func<double, double> weight = t => Math.Exp(t - Math.Exp(t));
            // Split near the mass of log(Exp(1)); a single wide Simpson panel can miss it.
            double mean = Integrate(t => value(t) * weight(t), -100, 0) + Integrate(t => value(t) * weight(t), 0, 4);
            if (!variance) return mean / scale / rate;
            Func<double, double> centered = t => { double d = value(t) - mean; return d * d * weight(t); };
            return (Integrate(centered, -100, 0) + Integrate(centered, 0, 4)) / scale / scale / rate / rate;
        }
    }
}
