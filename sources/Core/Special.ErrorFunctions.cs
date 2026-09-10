using System;
using System.Numerics;

namespace UMapx.Core
{
    public static partial class Special
    {
        private static readonly double[] GaussNodes =
        {
            0.0950125098376374402, 0.281603550779258913, 0.458016777657227386, 0.617876244402643748,
            0.755404408355003034, 0.865631202387831744, 0.944575023073232576, 0.989400934991649933
        };
        private static readonly double[] GaussWeights =
        {
            0.189450610455068496, 0.182603415044923589, 0.169156519395002538, 0.149595988816576733,
            0.124628971255533872, 0.0951585116824927848, 0.0622535239386478929, 0.0271524594117540949
        };

        private static Complex IntegrateGauss(Func<double, Complex> f, double end, int panels)
        {
            Complex sum = 0;
            double half = end / (2 * panels);
            for (int j = 0; j < panels; j++)
            {
                double mid = (2 * j + 1) * half;
                for (int k = 0; k < GaussNodes.Length; k++)
                    sum += GaussWeights[k] * (f(mid - half * GaussNodes[k]) + f(mid + half * GaussNodes[k]));
            }
            return sum * half;
        }

        // Fourier-Laplace representation in the upper half-plane. Reflection is
        // performed before evaluation, so the integral always has a decaying kernel.
        private static Complex FaddeevaValue(Complex z)
        {
            if (!IsFinite(z)) return ComplexNaN;
            if (z.Imaginary < 0) return 2 * Complex.Exp(-z * z) - FaddeevaValue(-z);
            if (z.Magnitude >= 12)
            {
                Complex term = 1, sum = 1;
                double previous = double.PositiveInfinity;
                for (int n = 1; n < 200; n++)
                {
                    term *= (n - 0.5) / (z * z);
                    if (term.Magnitude > previous) break;
                    sum += term;
                    if (term.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
                    previous = term.Magnitude;
                }
                Complex result = Complex.ImaginaryOne * sum / (SqrtPi * z);
                // The exponentially small real part on the real axis is not in
                // the algebraic expansion, but is still representable near the switch.
                return z.Imaginary == 0 ? new Complex(Math.Exp(-z.Real * z.Real), result.Imaginary) : result;
            }
            return IntegrateGauss(t => Complex.Exp(new Complex(-t * t / 4 - z.Imaginary * t, z.Real * t)),
                12, Math.Max(12, (int)Math.Ceiling(2 * z.Magnitude))) / SqrtPi;
        }

        private static double ErfcValue(double x)
        {
            if (double.IsNaN(x)) return double.NaN;
            if (x < 0) return 2 - ErfcValue(-x);
            if (double.IsPositiveInfinity(x)) return 0;
            return IncompleteGamma(0.5, x * x, true, true);
        }

        private static double ErfValue(double x)
        {
            if (double.IsNaN(x)) return double.NaN;
            if (double.IsInfinity(x)) return Math.Sign(x);
            if (Math.Abs(x) > 0.5) return x < 0 ? ErfcValue(-x) - 1 : 1 - ErfcValue(x);
            double term = x, sum = x;
            for (int n = 1; n < 100; n++)
            {
                term *= -x * x / n;
                double add = term / (2 * n + 1);
                sum += add;
                if (Math.Abs(add) <= SpecialEpsilon * Math.Abs(sum)) break;
            }
            return 2 * sum / SqrtPi;
        }

        private static Complex ErfcValue(Complex z)
        {
            if (z.Imaginary == 0) return new Complex(ErfcValue(z.Real), 0);
            if (z.Real < 0) return 2 - ErfcValue(-z);
            return Complex.Exp(-z * z) * FaddeevaValue(Complex.ImaginaryOne * z);
        }

        private static Complex ErfValue(Complex z)
        {
            if (z.Imaginary == 0) return new Complex(ErfValue(z.Real), 0);
            if (!IsFinite(z)) return ComplexNaN;
            if (z.Real < 0) return -ErfValue(-z);
            if (z.Magnitude > 0.5) return 1 - ErfcValue(z);
            Complex term = z, sum = z;
            for (int n = 1; n < 100; n++)
            {
                term *= -z * z / n;
                Complex add = term / (2 * n + 1);
                sum += add;
                if (add.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
            }
            return 2 * sum / SqrtPi;
        }

        private static double InverseErfc(double p)
        {
            if (double.IsNaN(p) || p < 0 || p > 2) return double.NaN;
            if (p == 0) return double.PositiveInfinity;
            if (p == 2) return double.NegativeInfinity;
            if (p == 1) return 0;
            if (p > 1) return -InverseErfc(2 - p);
            double lower = 0, upper = 28;
            // Bracketing also resolves very small probabilities without forming 1-p.
            for (int j = 0; j < 60; j++)
            {
                double mid = (lower + upper) / 2;
                if (ErfcValue(mid) > p) lower = mid; else upper = mid;
            }
            return (lower + upper) / 2;
        }

        private static double InverseErf(double x)
        {
            if (double.IsNaN(x) || Math.Abs(x) > 1) return double.NaN;
            if (Math.Abs(x) < 1e-4) return SqrtPi / 2 * x * (1 + Math.PI * x * x / 12);
            return x < 0 ? -InverseErfc(1 + x) : InverseErfc(1 - x);
        }

        private static Complex InverseErf(Complex z)
        {
            if (z.Imaginary == 0 && Math.Abs(z.Real) <= 1) return new Complex(InverseErf(z.Real), 0);
            if (!IsFinite(z)) return ComplexNaN;
            if (z.Real < 0) return -InverseErf(-z);
            Complex log = Complex.Log(1 - z * z);
            Complex t = 2 / (Math.PI * 0.147) + log / 2;
            Complex w = z.Magnitude < 0.5 ? SqrtPi * z / 2 : Complex.Sqrt(Complex.Sqrt(t * t - log / 0.147) - t);
            for (int n = 0; n < 50; n++)
            {
                Complex correction = (ErfValue(w) - z) / (2 / SqrtPi * Complex.Exp(-w * w));
                correction /= 1 + w * correction;
                w -= correction;
                if (correction.Magnitude <= SpecialEpsilon * (1 + w.Magnitude)) return w;
            }
            return ComplexNaN;
        }

        private static Complex DawsonValue(Complex z)
        {
            if (!IsFinite(z)) return ComplexNaN;
            if (z.Imaginary < 0) return -DawsonValue(-z);
            if (z.Magnitude < 0.5)
            {
                Complex term = z, sum = z;
                for (int n = 1; n < 100; n++)
                {
                    term *= -2 * z * z / (2 * n + 1);
                    sum += term;
                    if (term.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
                }
                return sum;
            }
            return -Complex.ImaginaryOne * SqrtPi / 2 * (FaddeevaValue(z) - Complex.Exp(-z * z));
        }

        private static Complex FresnelValue(Complex z, bool sine)
        {
            if (z.Magnitude < 1)
            {
                Complex term = sine ? z * z * z / 3 : z, sum = term;
                for (int k = 0; k < 100; k++)
                {
                    term *= sine ? -(4 * k + 3) * z * z * z * z / ((2.0 * k + 3) * (2 * k + 2) * (4 * k + 7))
                        : -(4 * k + 1) * z * z * z * z / ((2.0 * k + 2) * (2 * k + 1) * (4 * k + 5));
                    sum += term;
                    if (term.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
                }
                return sum;
            }
            Complex rotation = new Complex(Math.Sqrt(0.5), Math.Sqrt(0.5));
            Complex positive = SqrtPi / 2 * rotation * ErfValue(Complex.Conjugate(rotation) * z);
            if (z.Imaginary == 0) return sine ? positive.Imaginary : positive.Real;
            Complex negative = SqrtPi / 2 * Complex.Conjugate(rotation) * ErfValue(rotation * z);
            return sine ? (positive - negative) / (2 * Complex.ImaginaryOne) : (positive + negative) / 2;
        }

        private static Complex GeneralizedErf(Complex z, int n)
        {
            if (n < 0 || !IsFinite(z)) return ComplexNaN;
            if (z == Complex.Zero) return Complex.Zero;
            if (n == 0) return z / (Math.E * SqrtPi);
            if (n == 2) return ErfValue(z);
            Complex power = Complex.Pow(z, n);
            if (power.Magnitude < 8)
            {
                Complex term = 1, sum = 1;
                for (int k = 1; k < 1000; k++)
                {
                    term *= -power / k;
                    Complex add = term / ((double)n * k + 1);
                    sum += add;
                    if (add.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
                }
                return GammaValue((double)n + 1) / SqrtPi * z * sum;
            }
            // The power correction preserves the entire continuation in z when
            // z^n crosses the principal cut of the incomplete gamma function.
            Complex correction = z / Complex.Pow(power, 1.0 / n);
            return GammaValue((double)n) / SqrtPi * correction * IncompleteGamma((Complex)(1.0 / n), power, false, false);
        }

        private static Complex LambertValue(Complex z, int branch)
        {
            if (!IsFinite(z)) return ComplexNaN;
            if (z == Complex.Zero) return branch == 0 ? Complex.Zero : ComplexNaN;
            Complex w;
            Complex p = Complex.Sqrt(2 * (Math.E * z + 1));
            bool nearBranch = (Math.E * z + 1).Magnitude < 0.7;
            if (nearBranch && (branch == 0 || (branch == -1 && z.Imaginary >= 0) || (branch == 1 && z.Imaginary < 0)))
            {
                if (branch != 0) p = -p;
                w = -1 + p - p * p / 3 + 11 * p * p * p / 72;
            }
            else if (branch == 0 && z.Magnitude < 3) w = Complex.Log(1 + z);
            else if (branch == -1 && z.Imaginary == 0 && z.Real > -1 / Math.E && z.Real < 0)
            {
                double log = Math.Log(-z.Real);
                w = log - Math.Log(-log);
            }
            else
            {
                Complex log = Complex.Log(z) + new Complex(0, 2 * Math.PI * branch);
                w = log - Complex.Log(log) + Complex.Log(log) / log;
            }
            for (int n = 0; n < 100; n++)
            {
                Complex f = w - z * Complex.Exp(-w);
                if (f.Magnitude <= SpecialEpsilon * w.Magnitude) return w;
                Complex step = f / (w + 1 - (w + 2) * f / (2 * w + 2));
                w -= step;
                if (step.Magnitude <= SpecialEpsilon * (1 + w.Magnitude)) return w;
                if (!IsFinite(w)) return ComplexNaN;
            }
            return ComplexNaN;
        }
    }
}
