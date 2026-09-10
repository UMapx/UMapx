using System;
using System.Numerics;

namespace UMapx.Core
{
    public static partial class Special
    {
        // Keep intermediate arithmetic in double precision. Public single-precision
        // overloads round only once, after the tail or logarithmic ratio is evaluated.
        private const double SpecialEpsilon = 2e-15;
        private const double SqrtPi = 1.7724538509055160273;
        private static readonly Complex ComplexNaN = new Complex(double.NaN, double.NaN);
        private static readonly double[] GammaCoefficients =
        {
            0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
        };

        private static bool IsPole(double x) => x <= 0 && x == Math.Floor(x);
        private static bool IsPole(Complex x) => x.Imaginary == 0 && IsPole(x.Real);
        private static bool IsFinite(Complex x) => !double.IsNaN(x.Real) && !double.IsInfinity(x.Real)
            && !double.IsNaN(x.Imaginary) && !double.IsInfinity(x.Imaginary);

        private static double GammaLog(double x)
        {
            if (double.IsNaN(x) || IsPole(x)) return double.NaN;
            if (double.IsPositiveInfinity(x)) return x;
            if (x < 0.5)
                return Math.Log(Math.PI / Math.Abs(Math.Sin(Math.PI * (x % 2)))) - GammaLog(1 - x);
            double z = x - 1, sum = GammaCoefficients[0];
            for (int k = 1; k < GammaCoefficients.Length; k++) sum += GammaCoefficients[k] / (z + k);
            double t = z + 7.5;
            return 0.91893853320467274178 + Math.Log(sum) + (z + 0.5) * Math.Log(t) - t;
        }

        private static double GammaSign(double x) => x > 0 ? 1 : Math.Sign(Math.Sin(Math.PI * (x % 2)));

        private static double GammaValue(double x)
        {
            double log = GammaLog(x);
            return double.IsNaN(log) ? double.NaN : GammaSign(x) * Math.Exp(log);
        }

        // Analytic log-gamma on the plane cut along the negative real axis.
        // Recurrence retains the winding of log Gamma instead of taking log(Gamma).
        private static Complex GammaLog(Complex x)
        {
            if (!IsFinite(x) || IsPole(x)) return ComplexNaN;
            if (x.Real < -128)
            {
                double sign = x.Imaginary < 0 ? -1 : 1;
                Complex reflected = Math.Log(Math.PI) - LogSinPi(x) - GammaLog(1 - x);
                double turns = Math.Floor((1 - x.Real) / 2);
                return reflected - new Complex(0, sign * 2 * Math.PI * turns);
            }
            Complex correction = Complex.Zero;
            while (x.Real < 0.5)
            {
                correction -= Complex.Log(x);
                x += 1;
            }
            Complex z = x - 1, sum = GammaCoefficients[0];
            for (int k = 1; k < GammaCoefficients.Length; k++) sum += GammaCoefficients[k] / (z + k);
            Complex t = z + 7.5;
            return correction + 0.91893853320467274178 + Complex.Log(sum) + (z + 0.5) * Complex.Log(t) - t;
        }

        private static Complex GammaValue(Complex x) => Complex.Exp(GammaLog(x));

        private static double GammaSeries(double s, double x)
        {
            double term = 1 / s, sum = term;
            for (int n = 1; n <= 100000; n++)
            {
                term *= x / (s + n);
                sum += term;
                if (Math.Abs(term) <= Math.Abs(sum) * SpecialEpsilon) return sum;
            }
            return double.NaN;
        }

        private static Complex GammaSeries(Complex s, Complex x)
        {
            Complex term = 1 / s, sum = term;
            for (int n = 1; n <= 100000; n++)
            {
                term *= x / (s + n);
                sum += term;
                if (term.Magnitude <= sum.Magnitude * SpecialEpsilon) return sum;
                if (!IsFinite(sum)) return ComplexNaN;
            }
            return ComplexNaN;
        }

        private static double GammaFraction(double s, double x)
        {
            const double tiny = 1e-290;
            double b = x + 1 - s, c = 1 / tiny, d = 1 / (Math.Abs(b) < tiny ? tiny : b), h = d;
            for (int i = 1; i <= 100000; i++)
            {
                double a = i * (s - i);
                b += 2;
                d = b + a * d;
                c = b + a / c;
                if (Math.Abs(d) < tiny) d = tiny;
                if (Math.Abs(c) < tiny) c = tiny;
                d = 1 / d;
                double delta = c * d;
                h *= delta;
                if (Math.Abs(delta - 1) <= SpecialEpsilon) return h;
            }
            return double.NaN;
        }

        private static Complex GammaFraction(Complex s, Complex x)
        {
            const double tiny = 1e-290;
            Complex b = x + 1 - s, c = 1 / tiny, d = 1 / (b.Magnitude < tiny ? tiny : b), h = d;
            for (int i = 1; i <= 100000; i++)
            {
                Complex a = i * (s - i);
                b += 2;
                d = b + a * d;
                c = b + a / c;
                if (d.Magnitude < tiny) d = tiny;
                if (c.Magnitude < tiny) c = tiny;
                d = 1 / d;
                Complex delta = c * d;
                h *= delta;
                if ((delta - 1).Magnitude <= SpecialEpsilon) return h;
                if (!IsFinite(h)) return ComplexNaN;
            }
            return ComplexNaN;
        }

        private static double IncompleteGamma(double s, double x, bool upper, bool regularized)
        {
            if (double.IsNaN(s) || double.IsNaN(x) || double.IsInfinity(s) || IsPole(s) || x < 0) return double.NaN;
            if (x == 0) return upper ? (regularized ? 1 : GammaValue(s)) : (s > 0 ? 0 : double.NaN);
            if (double.IsPositiveInfinity(x)) return upper ? 0 : (regularized ? 1 : GammaValue(s));
            if (s > 0 && s < 0.01 && x <= 1)
            {
                // Evaluate Q directly when P is almost one. The usual subtraction
                // loses the entire tail for shapes close to zero.
                double logGammaOnePlus = -EulerGamma * s;
                double power = s * s;
                for (int k = 2; k <= 10; k++)
                {
                    logGammaOnePlus += (k % 2 == 0 ? 1 : -1) * ZetaValue((Complex)k).Real * power / k;
                    power *= s;
                }
                double log = s * Math.Log(x) - logGammaOnePlus;
                double term = -x, sum = s * term / (s + 1);
                for (int k = 2; k < 1000; k++)
                {
                    term *= -x / k;
                    double add = s * term / (s + k);
                    sum += add;
                    if (Math.Abs(add) <= SpecialEpsilon * Math.Abs(sum)) break;
                }
                double q = -Expm1(log) - Math.Exp(log) * sum;
                double value = upper ? q : 1 - q;
                return regularized ? value : value * GammaValue(s);
            }
            double scale = s * Math.Log(x) - x;
            double total = regularized ? 1 : GammaValue(s);
            double sign = 1;
            if (regularized) { scale -= GammaLog(s); sign = GammaSign(s); }
            if (x < s + 1)
            {
                double series = GammaSeries(s, x);
                if (double.IsNaN(series)) return double.NaN;
                double lower = sign * Math.Sign(series) * Math.Exp(scale + Math.Log(Math.Abs(series)));
                return upper ? total - lower : lower;
            }
            double tail = sign * Math.Exp(scale) * GammaFraction(s, x);
            return upper ? tail : total - tail;
        }

        private static Complex IncompleteGamma(Complex s, Complex x, bool upper, bool regularized)
        {
            if (s.Imaginary == 0 && x.Imaginary == 0 && x.Real >= 0)
                return new Complex(IncompleteGamma(s.Real, x.Real, upper, regularized), 0);
            if (!IsFinite(s) || !IsFinite(x) || IsPole(s)) return ComplexNaN;
            if (x == Complex.Zero) return upper ? (regularized ? Complex.One : GammaValue(s)) : (s.Real > 0 ? Complex.Zero : ComplexNaN);
            Complex scale = s * Complex.Log(x) - x - (regularized ? GammaLog(s) : Complex.Zero);
            Complex total = regularized ? Complex.One : GammaValue(s);
            if (x.Real < 0 || x.Magnitude < s.Magnitude + 1)
            {
                Complex lower = Complex.Exp(scale) * GammaSeries(s, x);
                return upper ? total - lower : lower;
            }
            Complex tail = Complex.Exp(scale) * GammaFraction(s, x);
            return upper ? tail : total - tail;
        }

        private static double BetaLog(double a, double b) => GammaLog(a) + GammaLog(b) - GammaLog(a + b);
        private static Complex BetaLog(Complex a, Complex b) => GammaLog(a) + GammaLog(b) - GammaLog(a + b);
        private static double BetaValue(double a, double b)
        {
            double log = BetaLog(a, b);
            return double.IsNaN(log) ? double.NaN : GammaSign(a) * GammaSign(b) * GammaSign(a + b) * Math.Exp(log);
        }

        private static double BetaFraction(double a, double b, double x)
        {
            const double tiny = 1e-290;
            double c = 1, d = 1 - (a + b) * x / (a + 1);
            if (Math.Abs(d) < tiny) d = tiny;
            d = 1 / d;
            double h = d;
            for (int m = 1; m <= 100000; m++)
            {
                double aa = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m));
                d = 1 + aa * d; c = 1 + aa / c;
                if (Math.Abs(d) < tiny) d = tiny;
                if (Math.Abs(c) < tiny) c = tiny;
                d = 1 / d; h *= d * c;
                aa = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1));
                d = 1 + aa * d; c = 1 + aa / c;
                if (Math.Abs(d) < tiny) d = tiny;
                if (Math.Abs(c) < tiny) c = tiny;
                d = 1 / d;
                double delta = d * c;
                h *= delta;
                if (Math.Abs(delta - 1) <= SpecialEpsilon) return h;
            }
            return double.NaN;
        }

        private static double IncompleteBeta(double a, double b, double x, bool regularized)
        {
            if (!(a > 0) || !(b > 0) || double.IsInfinity(a + b) || double.IsNaN(x) || x < 0 || x > 1) return double.NaN;
            if (x == 0) return 0;
            double logBeta = BetaLog(a, b);
            double total = regularized ? 1 : Math.Exp(logBeta);
            if (x == 1) return total;
            bool complement = x > (a + 1) / (a + b + 2);
            if (complement) { double t = a; a = b; b = t; x = 1 - x; }
            double value = Math.Exp(a * Math.Log(x) + b * Math.Log(1 - x) - (regularized ? logBeta : 0)) * BetaFraction(a, b, x) / a;
            return complement ? total - value : value;
        }

        private static Complex ChebyshevValue(Complex x, int n, bool secondKind)
        {
            if (n == int.MinValue) return ComplexNaN;
            if (n < 0) return secondKind ? (n == -1 ? Complex.Zero : -ChebyshevValue(x, -n - 2, true)) : ChebyshevValue(x, -n, false);
            if (x == Complex.One) return secondKind ? (Complex)((double)n + 1) : Complex.One;
            if (x == -Complex.One) return ((n & 1) == 0 ? 1 : -1) * (secondKind ? (Complex)((double)n + 1) : Complex.One);
            if (n > 100000)
            {
                Complex angle = Complex.Acos(x);
                return secondKind ? Complex.Sin(((double)n + 1) * angle) / Complex.Sin(angle) : Complex.Cos(n * angle);
            }
            Complex previous = 1, current = secondKind ? 2 * x : x;
            if (n == 0) return previous;
            for (int k = 1; k < n; k++) { Complex next = 2 * x * current - previous; previous = current; current = next; }
            return current;
        }

        private static Complex FactorialProduct(Complex n, Complex k, bool rising)
        {
            if (k == Complex.Zero) return Complex.One;
            if (k.Imaginary == 0 && k.Real >= 0 && k.Real <= 10000 && k.Real == Math.Floor(k.Real))
            {
                Complex product = 1;
                for (int j = 0; j < (int)k.Real; j++)
                {
                    Complex factor = rising ? n + j : n - j;
                    if (factor == Complex.Zero) return Complex.Zero;
                    product *= factor;
                }
                return product;
            }
            return rising ? Complex.Exp(GammaLog(n + k) - GammaLog(n)) : Complex.Exp(GammaLog(n + 1) - GammaLog(n - k + 1));
        }

        private static Complex BinomialValue(Complex n, Complex k)
        {
            if (k.Imaginary == 0 && k.Real < 0) return Complex.Zero;
            if (k == Complex.Zero) return Complex.One;
            if (n.Imaginary == 0 && k.Imaginary == 0 && n.Real >= 0 && n.Real == Math.Floor(n.Real) && k.Real == Math.Floor(k.Real))
            {
                if (k.Real > n.Real) return Complex.Zero;
                k = Math.Min(k.Real, n.Real - k.Real);
            }
            if (k.Imaginary == 0 && k.Real >= 0 && k.Real <= 10000 && k.Real == Math.Floor(k.Real))
            {
                Complex product = 1;
                for (int j = 1; j <= (int)k.Real; j++) product *= (n - j + 1) / j;
                return product;
            }
            return Complex.Exp(GammaLog(n + 1) - GammaLog(k + 1) - GammaLog(n - k + 1));
        }

        private static Complex Hypergeometric2F1(Complex a, Complex b, Complex c, Complex z)
        {
            if (z == Complex.Zero || a == Complex.Zero || b == Complex.Zero) return Complex.One;
            if (IsPole(c)) return ComplexNaN;
            bool terminating = IsPole(a) || IsPole(b);
            if (terminating && (1 - z).Magnitude < 0.5 && (c - a - b).Real > 0 && (IsPole(c - a) || IsPole(c - b)))
                return Complex.Pow(1 - z, c - a - b) * Hypergeometric2F1(c - a, c - b, c, z);
            if (!terminating && z.Real < 0)
                return Complex.Pow(1 - z, -a) * Hypergeometric2F1(a, c - b, c, z / (z - 1));
            if (!terminating && z == Complex.One)
                return (c - a - b).Real > 0 ? Complex.Exp(GammaLog(c) + GammaLog(c - a - b) - GammaLog(c - a) - GammaLog(c - b)) : ComplexNaN;
            if (!terminating && z.Magnitude > 0.8) return HypergeometricContinue(a, b, c, z, 2);
            Complex term = 1, sum = 1;
            for (int n = 0; n < 100000; n++)
            {
                term *= ((a + n) / (c + n)) * ((b + n) / (n + 1)) * z;
                sum += term;
                if (term.Magnitude <= SpecialEpsilon * sum.Magnitude || term == Complex.Zero) return sum;
                if (!IsFinite(sum)) return ComplexNaN;
            }
            return ComplexNaN;
        }
    }
}
