using System;
using System.Numerics;

namespace UMapx.Core
{
    public static partial class Special
    {
        private static Complex OrthogonalPolynomial(Complex x, Complex a, int n, int family)
        {
            if (n < 0) return ComplexNaN;
            Complex previous = 1;
            if (n == 0) return previous;
            Complex current = family == 0 ? 1 + a - x : family == 1 ? x : family == 2 ? 2 * x : 2 * a * x;
            for (int k = 2; k <= n; k++)
            {
                Complex next;
                switch (family)
                {
                    case 0: next = ((2 * k - 1 + a - x) * current - (k - 1 + a) * previous) / k; break;
                    case 1: next = ((2 * k - 1) * x * current - (k - 1) * previous) / k; break;
                    case 2: next = 2 * (x * current - (k - 1) * previous); break;
                    default: next = (2 * x * (k + a - 1) * current - (k + 2 * a - 2) * previous) / k; break;
                }
                previous = current; current = next;
            }
            return current;
        }

        private static Complex ErlangBlocking(Complex y, int n)
        {
            if (n < 0 || !IsFinite(y)) return ComplexNaN;
            Complex blocking = 1;
            for (int k = 1; k <= n; k++) blocking = y * blocking / (k + y * blocking);
            return blocking;
        }

        private static Complex LogSinPi(Complex z)
        {
            double real = z.Real % 2;
            if (Math.Abs(z.Imaginary) < 20) return Complex.Log(Complex.Sin(Math.PI * new Complex(real, z.Imaginary)));
            double phase = (z.Imaginary > 0 ? 1 : -1) * (Math.PI / 2 - Math.PI * real);
            return new Complex(Math.PI * Math.Abs(z.Imaginary) - Math.Log(2), Math.Atan2(Math.Sin(phase), Math.Cos(phase)));
        }

        private static Complex Polygamma(Complex z, bool derivative)
        {
            if (!IsFinite(z) || IsPole(z)) return ComplexNaN;
            if (z.Real < 0.5)
            {
                Complex reduced = new Complex(z.Real % 1, z.Imaginary);
                if (derivative)
                    return Math.PI * Math.PI * Complex.Exp(-2 * LogSinPi(z)) - Polygamma(1 - z, true);
                Complex cotangent = Math.Abs(z.Imaginary) > 100 ? new Complex(0, z.Imaginary > 0 ? -1 : 1)
                    : Complex.Cos(Math.PI * reduced) / Complex.Sin(Math.PI * reduced);
                return Polygamma(1 - z, false) - Math.PI * cotangent;
            }
            Complex accumulated = 0;
            while (z.Real < 12)
            {
                accumulated += derivative ? 1 / (z * z) : -1 / z;
                z += 1;
            }
            Complex inverse = 1 / z, square = inverse * inverse;
            Complex result = derivative ? inverse + square / 2 : Complex.Log(z) - inverse / 2;
            double[] bernoulli = { 1.0 / 6, -1.0 / 30, 1.0 / 42, -1.0 / 30, 5.0 / 66, -691.0 / 2730, 7.0 / 6, -3617.0 / 510 };
            Complex power = square;
            for (int k = 0; k < bernoulli.Length; k++)
            {
                result += derivative ? bernoulli[k] * power * inverse : -bernoulli[k] * power / (2 * k + 2);
                power *= square;
            }
            return accumulated + result;
        }

        private static Complex ZetaValue(Complex s)
        {
            if (!IsFinite(s) || s == Complex.One) return ComplexNaN;
            if (s == Complex.Zero) return -0.5;
            if (s.Real < 0)
            {
                if (s.Imaginary == 0 && s.Real % 2 == 0) return Complex.Zero;
                return Complex.Exp(s * Math.Log(2) + (s - 1) * Math.Log(Math.PI) + LogSinPi(s / 2) + GammaLog(1 - s)) * ZetaValue(1 - s);
            }
            if (s.Real > 55) return Complex.One;
            int count = (int)Math.Min(100000, 24 + Math.Ceiling(s.Imaginary == 0 ? 0 : Math.Abs(s.Imaginary)));
            Complex sum = 0;
            for (int k = 1; k < count; k++) sum += Complex.Exp(-s * Math.Log(k));
            Complex power = Complex.Exp(-s * Math.Log(count));
            sum += power * (count / (s - 1) + 0.5);
            double[] coefficients = { 1.0 / 12, -1.0 / 720, 1.0 / 30240, -1.0 / 1209600,
                1.0 / 47900160, -691.0 / 1307674368000, 1.0 / 74724249600, -3617.0 / 10670622842880000 };
            Complex term = s * power / count;
            for (int k = 0; k < coefficients.Length; k++)
            {
                sum += coefficients[k] * term;
                term *= (s + 2 * k + 1) / count * ((s + 2 * k + 2) / count);
            }
            return sum;
        }

        // Continue a hypergeometric solution by Taylor steps inside disks that
        // exclude the differential equation's singularities. Values and first
        // derivatives are propagated together, without subtracting large gamma ratios.
        private static void HypergeometricStep(Complex a, Complex b, Complex c, Complex center,
            Complex step, int kind, ref Complex value, ref Complex derivative)
        {
            Complex previous = value, current = step * derivative;
            Complex sum = previous + current, slope = current;
            for (int n = 0; n < 100; n++)
            {
                Complex next;
                if (kind == 2)
                    next = ((a + n) * (b + n) * step * step * previous - (n + 1) * (c - (a + b + 1) * center + n * (1 - 2 * center)) * step * current)
                        / (center * (1 - center) * (n + 1) * (n + 2));
                else
                    next = ((kind == 1 ? a + n : Complex.One) * step * step * previous - (n + 1) * (b + n - (kind == 1 ? center : Complex.Zero)) * step * current)
                        / (center * (n + 1) * (n + 2));
                sum += next;
                slope += (n + 2) * next;
                previous = current; current = next;
                if (n > 6 && next.Magnitude * (n + 2) <= SpecialEpsilon * (sum.Magnitude + slope.Magnitude)) break;
            }
            value = sum; derivative = slope / step;
        }

        private static Complex HypergeometricContinue(Complex a, Complex b, Complex c, Complex target, int kind)
        {
            Complex waypoint = kind == 2 && target.Real > 1 && Math.Abs(target.Imaginary) < 0.5
                ? new Complex(0.5, target.Imaginary > 0 ? 0.5 : -0.5) : target;
            Complex center = waypoint / waypoint.Magnitude * 0.25;
            Complex term = 1, value = 1, derivative = 0;
            for (int n = 0; n < 1000; n++)
            {
                term *= kind == 2 ? (a + n) / (c + n) * ((b + n) / (n + 1)) * center
                    : (kind == 1 ? a + n : Complex.One) / (b + n) * center / (n + 1);
                value += term; derivative += (n + 1) * term / center;
                if (term.Magnitude <= SpecialEpsilon * value.Magnitude) break;
            }
            for (int stage = 0; stage < 2; stage++)
            {
                Complex end = stage == 0 ? waypoint : target;
                for (int iteration = 0; iteration < 10000 && center != end; iteration++)
                {
                    double radius = kind == 2 ? Math.Min(center.Magnitude, (1 - center).Magnitude) : center.Magnitude;
                    double length = Math.Min(4, 0.4 * radius);
                    Complex delta = end - center;
                    Complex step = delta.Magnitude <= length ? delta : delta / delta.Magnitude * length;
                    HypergeometricStep(a, b, c, center, step, kind, ref value, ref derivative);
                    center = step == delta ? end : center + step;
                    if (!IsFinite(value)) return value;
                }
                if (center != end) return ComplexNaN;
            }
            return value;
        }

        private static Complex Hypergeometric1F1(Complex a, Complex b, Complex z)
        {
            if (z == Complex.Zero) return Complex.One;
            bool absentA = double.IsNaN(a.Real), absentB = double.IsNaN(b.Real);
            if ((absentA && absentB) || a == b) return Complex.Exp(z);
            if (absentB) return Complex.Pow(1 - z, -a);
            if (IsPole(b)) return ComplexNaN;
            if (!absentA && a == Complex.Zero) return Complex.One;
            if (z.Magnitude > 8) return HypergeometricContinue(a, b, Complex.Zero, z, absentA ? 0 : 1);
            Complex term = 1, sum = 1;
            for (int n = 0; n < 10000; n++)
            {
                term *= (absentA ? Complex.One : a + n) / (b + n) * z / (n + 1);
                sum += term;
                if (term.Magnitude <= SpecialEpsilon * sum.Magnitude || term == Complex.Zero) return sum;
            }
            return ComplexNaN;
        }
    }
}
