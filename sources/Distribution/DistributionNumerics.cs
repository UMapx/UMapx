using System;
using UMapx.Core;

namespace UMapx.Distribution
{
    internal static class DistributionNumerics
    {
        internal const double EulerGamma = 0.57721566490153286061;
        internal const double LogSqrtTwoPi = 0.91893853320467274178;

        internal static double Log1p(double x)
        {
            double sum = 1 + x;
            if (sum == 1) return x;
            if (double.IsPositiveInfinity(sum)) return sum;
            return Math.Log(sum) * x / (sum - 1);
        }

        internal static double Softplus(double x) => x > 0 ? x + Log1p(Math.Exp(-x)) : Log1p(Math.Exp(x));

        internal static double LogNormalSurvival(double x)
        {
            if (x < 0) return Log1p(-0.5 * Special.DistributionErfc(-x / Math.Sqrt(2)));
            if (x < 26) return Math.Log(0.5 * Special.DistributionErfc(x / Math.Sqrt(2)));
            return -0.5 * x * x + LogNormalMillsRatio(x);
        }

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

        internal static double BetaPrimeExcess(double a, double b)
        {
            if (b <= 4) return double.PositiveInfinity;
            return 6 * (a * (a + b - 1) * (5 * b - 11) + (b - 1) * (b - 1) * (b - 2)) /
                (a * (a + b - 1) * (b - 3) * (b - 4));
        }

        internal static double PoissonLogMass(double lambda, double k)
        {
            return k == 0 ? -lambda : -lambda + k * Math.Log(lambda) - Special.DistributionLogGamma(k + 1);
        }

        // Adaptive integration on a finite transformed interval. The callers choose decaying tails.
        internal static double Integrate(Func<double, double> f, double a, double b, double tolerance = 1e-11)
        {
            double mid = (a + b) / 2, fa = f(a), fm = f(mid), fb = f(b);
            return Refine(f, a, b, fa, fm, fb, (b - a) * (fa + 4 * fm + fb) / 6, tolerance, 24);
        }

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
