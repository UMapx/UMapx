using System;
using System.Numerics;

namespace UMapx.Core
{
    public static partial class Special
    {
        private static double Expm1(double x)
        {
            if (Math.Abs(x) > 0.5) return Math.Exp(x) - 1;
            double term = x, sum = x;
            for (int n = 2; n < 100; n++)
            {
                term *= x / n; sum += term;
                if (Math.Abs(term) <= SpecialEpsilon * Math.Abs(sum)) break;
            }
            return sum;
        }

        private static Complex LogisticValue(Complex z)
        {
            if (z.Real >= 0) return 1 / (1 + Complex.Exp(-z));
            Complex exponential = Complex.Exp(z);
            return exponential / (1 + exponential);
        }

        private static Complex RademacherValue(Complex z, int n)
        {
            if (!IsFinite(z)) return ComplexNaN;
            if (n > 1023) return z.Imaginary == 0 ? Complex.Zero : new Complex(0, Math.Sign(z.Imaginary));
            Complex scaled = Math.Pow(2, n) * z;
            double real = scaled.Real % 2;
            if (z.Imaginary == 0)
                return real == Math.Truncate(real) ? Complex.Zero : (Complex)Math.Sign(Math.Sin(Math.PI * real));
            Complex sine;
            if (Math.Abs(scaled.Imaginary) > 20)
                sine = new Complex(Math.Sin(Math.PI * real), Math.Sign(scaled.Imaginary) * Math.Cos(Math.PI * real));
            else sine = Complex.Sin(Math.PI * new Complex(real, scaled.Imaginary));
            return sine / sine.Magnitude;
        }

        private static int FibonacciValue(int n, bool lucas)
        {
            int limit = lucas ? 44 : 46;
            if (n < -limit || n > limit) throw new ArgumentOutOfRangeException(nameof(n), "The result does not fit in Int32.");
            long previous = lucas ? 2 : 0, current = 1;
            int order = Math.Abs(n);
            for (int k = 0; k < order; k++) { long next = previous + current; previous = current; current = next; }
            if (n < 0 && ((order & 1) == (lucas ? 1 : 0))) previous = -previous;
            return (int)previous;
        }

        private static double NumberPolynomial(int n, double x, bool euler)
        {
            if (n < 0 || n > (euler ? 186 : 258)) return double.NaN;
            double variable = euler ? x - 0.5 : x;
            double sum = 0, binomial = 1, scale = 1;
            for (int k = 0; k <= n; k++)
            {
                sum = sum * variable + binomial * (euler ? Euler(k) * scale : Bernoulli(k));
                binomial *= (double)(n - k) / (k + 1);
                scale *= 0.5;
            }
            return sum;
        }
    }
}
