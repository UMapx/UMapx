using System;
using System.Numerics;

namespace UMapx.Core
{
    public static partial class Special
    {
        private const double EulerGamma = 0.57721566490153286061;

        private static Complex ImaginaryPower(int n)
        {
            switch ((n % 4 + 4) % 4)
            {
                case 0: return Complex.One;
                case 1: return Complex.ImaginaryOne;
                case 2: return -Complex.One;
                default: return -Complex.ImaginaryOne;
            }
        }

        private static Complex BesselSeries(Complex z, int n, bool modified)
        {
            Complex term = Complex.Exp(n * Complex.Log(z / 2) - GammaLog((double)n + 1));
            Complex sum = term, factor = (modified ? 1 : -1) * z * z / 4;
            for (int k = 1; k <= 10000; k++)
            {
                term *= factor / ((double)k * (n + k));
                sum += term;
                if (term.Magnitude <= SpecialEpsilon * sum.Magnitude) return sum;
                if (!IsFinite(sum)) return sum;
            }
            return ComplexNaN;
        }

        // Hankel expansion is evaluated only for orders zero and one. Higher
        // orders use recurrence; their first asymptotic correction need not be small.
        private static Complex BesselAsymptotic(Complex z, int n, bool secondKind)
        {
            Complex term = 1, even = 1, odd = 0;
            double previous = double.PositiveInfinity;
            for (int k = 1; k < 200; k++)
            {
                term *= (4.0 * n * n - (2.0 * k - 1) * (2.0 * k - 1)) / (8 * k * z);
                if (term.Magnitude > previous) break;
                if ((k & 1) == 0) even += ((k / 2 & 1) == 0 ? 1 : -1) * term;
                else odd += (((k - 1) / 2 & 1) == 0 ? 1 : -1) * term;
                if (term.Magnitude <= SpecialEpsilon * (even.Magnitude + odd.Magnitude)) break;
                previous = term.Magnitude;
            }
            // Separate the fixed phase shift from z to avoid losing it for large z.
            double phase = n * Math.PI / 2 + Math.PI / 4;
            Complex cosine = Complex.Cos(z) * Math.Cos(phase) + Complex.Sin(z) * Math.Sin(phase);
            Complex sine = Complex.Sin(z) * Math.Cos(phase) - Complex.Cos(z) * Math.Sin(phase);
            return Complex.Sqrt(2 / (Math.PI * z)) * (secondKind ? sine * even + cosine * odd : cosine * even - sine * odd);
        }

        private static Complex BesselJ(Complex z, int order)
        {
            if (!IsFinite(z) || order == int.MinValue) return ComplexNaN;
            int n = Math.Abs(order);
            if (n > 100000) return ComplexNaN;
            double sign = order < 0 && (n & 1) != 0 ? -1 : 1;
            if (z == Complex.Zero) return n == 0 ? Complex.One : Complex.Zero;
            if (z.Real < 0) return ((n & 1) == 0 ? sign : -sign) * BesselJ(-z, n);
            if (z.Magnitude <= 12) return sign * BesselSeries(z, n, false);
            Complex j0 = BesselAsymptotic(z, 0, false), j1 = BesselAsymptotic(z, 1, false);
            if (n == 0) return j0;
            if (n == 1) return sign * j1;
            if (n <= (z.Imaginary == 0 ? z.Magnitude : z.Magnitude * 0.5))
            {
                for (int k = 1; k < n; k++) { Complex next = 2 * k / z * j1 - j0; j0 = j1; j1 = next; }
                return sign * j1;
            }
            // Miller recurrence selects the minimal solution J_n when upward
            // recurrence would amplify contamination by the dominant solution Y_n.
            int start = n + (int)Math.Min(100000, Math.Ceiling(z.Magnitude)) + 40;
            Complex current = 1, nextValue = 0, result = 0;
            for (int k = start; k >= 1; k--)
            {
                Complex previousValue = 2 * k / z * current - nextValue;
                nextValue = current;
                current = previousValue;
                if (k - 1 == n) result = current;
                if (current.Magnitude > 1e150)
                {
                    current *= 1e-150; nextValue *= 1e-150; result *= 1e-150;
                }
            }
            return sign * result * (j0.Magnitude >= j1.Magnitude ? j0 / current : j1 / nextValue);
        }

        private static Complex BesselYBase(Complex z, int n)
        {
            if (z.Magnitude > 12) return BesselAsymptotic(z, n, true);
            Complex term = n == 0 ? Complex.One : z / 2;
            double hk = 0, hnk = n;
            Complex sum = (hk + hnk - 2 * EulerGamma) * term;
            for (int k = 1; k < 10000; k++)
            {
                term *= -z * z / (4.0 * k * (n + k));
                hk += 1.0 / k; hnk += 1.0 / (n + k);
                Complex add = (hk + hnk - 2 * EulerGamma) * term;
                sum += add;
                if (add.Magnitude <= SpecialEpsilon * (1 + sum.Magnitude)) break;
            }
            return 2 / Math.PI * Complex.Log(z / 2) * BesselSeries(z, n, false) - sum / Math.PI - (n == 0 ? Complex.Zero : 2 / (Math.PI * z));
        }

        private static Complex BesselY(Complex z, int order)
        {
            if (!IsFinite(z) || order == int.MinValue) return ComplexNaN;
            int n = Math.Abs(order);
            if (n > 100000) return ComplexNaN;
            double sign = order < 0 && (n & 1) != 0 ? -1 : 1;
            if (z == Complex.Zero) return new Complex(sign * double.NegativeInfinity, 0);
            if (z.Real < 0)
                return sign * ((n & 1) == 0 ? 1 : -1) * (BesselY(-z, n) + (z.Imaginary < 0 ? -2 : 2) * Complex.ImaginaryOne * BesselJ(-z, n));
            if (z.Imaginary < 0) return sign * Complex.Conjugate(BesselY(Complex.Conjugate(z), n));
            if (z.Imaginary > 0)
                return sign * (Complex.ImaginaryOne * BesselJ(z, n) - 2 / Math.PI * ImaginaryPower(-n) * BesselK(-Complex.ImaginaryOne * z, n));
            Complex previous = BesselYBase(z, 0);
            if (n == 0) return previous;
            Complex current = BesselYBase(z, 1);
            for (int k = 1; k < n; k++) { Complex next = 2 * k / z * current - previous; previous = current; current = next; }
            return sign * current;
        }

        private static Complex BesselI(Complex z, int n)
        {
            if (n == int.MinValue) return ComplexNaN;
            n = Math.Abs(n);
            return ImaginaryPower(-n) * BesselJ(Complex.ImaginaryOne * z, n);
        }

        private static Complex BesselKBase(Complex z, int n)
        {
            if (z.Magnitude >= 9)
            {
                Complex term = 1, sum = 1;
                double previous = double.PositiveInfinity;
                for (int k = 1; k < 200; k++)
                {
                    term *= (4.0 * n * n - (2.0 * k - 1) * (2.0 * k - 1)) / (8 * k * z);
                    if (term.Magnitude > previous) break;
                    sum += term;
                    if (term.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
                    previous = term.Magnitude;
                }
                return Complex.Sqrt(Math.PI / (2 * z)) * Complex.Exp(-z) * sum;
            }
            Complex t = n == 0 ? Complex.One : z / 2;
            double hk = 0, hnk = n;
            Complex series = (hk + hnk - 2 * EulerGamma) * t;
            for (int k = 1; k < 10000; k++)
            {
                t *= z * z / (4.0 * k * (n + k));
                hk += 1.0 / k; hnk += 1.0 / (n + k);
                Complex add = (hk + hnk - 2 * EulerGamma) * t;
                series += add;
                if (add.Magnitude <= SpecialEpsilon * (1 + series.Magnitude)) break;
            }
            Complex logarithm = Complex.Log(z / 2) * BesselSeries(z, n, true);
            return n == 0 ? -logarithm + series / 2 : 1 / z + logarithm - series / 2;
        }

        private static Complex BesselK(Complex z, int order)
        {
            if (!IsFinite(z) || order == int.MinValue) return ComplexNaN;
            int n = Math.Abs(order);
            if (n > 100000) return ComplexNaN;
            if (z == Complex.Zero) return new Complex(double.PositiveInfinity, 0);
            if (z.Real < 0)
                return ((n & 1) == 0 ? 1 : -1) * BesselK(-z, n) - (z.Imaginary < 0 ? -1 : 1) * Complex.ImaginaryOne * Math.PI * BesselI(-z, n);
            Complex previous = BesselKBase(z, 0);
            if (n == 0) return previous;
            Complex current = BesselKBase(z, 1);
            for (int k = 1; k < n; k++) { Complex next = 2 * k / z * current + previous; previous = current; current = next; }
            return current;
        }

        private static Complex StruveValue(Complex z, int n, bool modified)
        {
            if (!IsFinite(z) || n < -100000 || n > 100000) return ComplexNaN;
            if (z == Complex.Zero) return n >= 0 ? Complex.Zero : n == -1 ? (Complex)(2 / Math.PI) : ComplexNaN;
            if (n < 0)
            {
                Complex previous = StruveValue(z, 0, modified);
                Complex current = 2 / Math.PI + (modified ? 1 : -1) * StruveValue(z, 1, modified);
                for (int k = -1; k > n; k--)
                {
                    Complex source = Complex.Pow(z / 2, k) / (SqrtPi * GammaValue(k + 1.5));
                    Complex next = 2 * k / z * current + (modified ? 1 : -1) * previous + source;
                    previous = current; current = next;
                }
                return current;
            }
            if (z.Magnitude < 4)
            {
                Complex term = Complex.Exp((n + 1) * Complex.Log(z / 2) - GammaLog(1.5) - GammaLog(n + 1.5));
                Complex sum = term;
                for (int k = 0; k < 1000; k++)
                {
                    term *= (modified ? 1 : -1) * z * z / (4 * (k + 1.5) * (k + n + 1.5));
                    sum += term;
                    if (term.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
                }
                return sum;
            }
            // The substitution t=cos(theta) removes the endpoint singularity in
            // the integral representation, including order zero (DLMF 11.5.1-2).
            int panels = Math.Max(4, (int)Math.Min(100000, Math.Ceiling(z.Magnitude / 4 + n / 4.0)));
            Complex integral = IntegrateGauss(theta => Math.Pow(Math.Sin(theta), 2.0 * n) *
                (modified ? Complex.Sinh(z * Math.Cos(theta)) : Complex.Sin(z * Math.Cos(theta))), Math.PI / 2, panels);
            return 2 / SqrtPi * Complex.Exp(n * Complex.Log(z / 2) - GammaLog(n + 0.5)) * integral;
        }
    }
}
