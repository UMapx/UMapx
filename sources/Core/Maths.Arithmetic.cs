using System;
using System.Collections.Generic;
using System.Numerics;

namespace UMapx.Core
{
    public static partial class Maths
    {
        private static Complex ComplexSinh(Complex32 value)
        {
            return new Complex(Math.Sinh(value.Real) * Math.Cos(value.Imag), Math.Cosh(value.Real) * Math.Sin(value.Imag));
        }

        private static Complex ComplexCosh(Complex32 value)
        {
            return new Complex(Math.Cosh(value.Real) * Math.Cos(value.Imag), Math.Sinh(value.Real) * Math.Sin(value.Imag));
        }

        private static double RealAsinh(double value)
        {
            double magnitude = Math.Abs(value);
            if (magnitude < 1e-4 || double.IsInfinity(value)) return value;
            double result = Math.Log(magnitude + Math.Sqrt(magnitude * magnitude + 1));
            return value < 0 ? -result : result;
        }

        private static double RealAcosh(double value)
        {
            if (value < 1) return double.NaN;
            return Math.Log(value + Math.Sqrt((value - 1) * (value + 1)));
        }

        private static double LogOnePlus(double value)
        {
            double sum = 1 + value;
            if (sum == 1) return value;
            return Math.Log(sum) * value / (sum - 1);
        }

        private static Complex ComplexLogOnePlus(Complex value)
        {
            double x = value.Real, y = value.Imaginary;
            return new Complex(0.5 * LogOnePlus(2 * x + x * x + y * y), Math.Atan2(y, 1 + x));
        }

        private static Complex PrincipalAtan(Complex value)
        {
            double x = value.Real, y = value.Imaginary;
            double denominator = x * x + (y - 1) * (y - 1);
            double ratio = 4 * y / denominator;
            double imaginary = Math.Abs(ratio) < 0.5 ? LogOnePlus(ratio) :
                Math.Log(x * x + (y + 1) * (y + 1)) - Math.Log(denominator);
            double real = x == 0 && Math.Abs(y) > 1 ? Math.Sign(y) * Math.PI / 2 :
                0.5 * Math.Atan2(2 * x, 1 - x * x - y * y);
            return new Complex(real, 0.25 * imaginary);
        }

        private static Complex32[] CubicWithSeparatedRoot(double a, double b, double c)
        {
            // When |a| dominates, depressing the cubic loses the two smaller roots.
            // Refine the isolated root and recover the remaining sum/product without cancellation.
            double root = -a;
            for (int i = 0; i < 4; i++)
                root -= (((root + a) * root + b) * root + c) / ((3 * root + 2 * a) * root + b);
            double product = -c / root;
            double sum = (b - product) / root;
            double discriminant = sum * sum - 4 * product;
            Complex first, second;
            if (discriminant < 0)
            {
                first = new Complex(sum / 2, Math.Sqrt(-discriminant) / 2);
                second = Complex.Conjugate(first);
            }
            else
            {
                double q = 0.5 * (sum + (sum < 0 ? -1 : 1) * Math.Sqrt(discriminant));
                first = q;
                second = q == 0 ? 0 : product / q;
            }
            return new[] { new Complex32((float)root, 0), (Complex32)first, (Complex32)second };
        }

        private static double RealAtanh(double value)
        {
            if (Math.Abs(value) < 1e-4) return value;
            return 0.5 * Math.Log((1 + value) / (1 - value));
        }

        private static double RealCubeRoot(double value)
        {
            double root = Math.Pow(Math.Abs(value), 1.0 / 3);
            return value < 0 ? -root : root;
        }

        private static ulong UnsignedMagnitude(long value)
        {
            return value < 0 ? (ulong)(-(value + 1)) + 1 : (ulong)value;
        }

        private static ulong UnsignedGcd(ulong a, ulong b)
        {
            while (b != 0)
            {
                ulong remainder = a % b;
                a = b;
                b = remainder;
            }
            return a;
        }

        private static BigInteger[] ExtendedGcd(BigInteger a, BigInteger b)
        {
            BigInteger oldR = a, r = b, oldS = 1, s = 0, oldT = 0, t = 1;
            while (r != 0)
            {
                BigInteger quotient = oldR / r;
                BigInteger nextR = oldR - quotient * r;
                BigInteger nextS = oldS - quotient * s;
                BigInteger nextT = oldT - quotient * t;
                oldR = r; r = nextR;
                oldS = s; s = nextS;
                oldT = t; t = nextT;
            }
            if (oldR < 0) { oldR = -oldR; oldS = -oldS; oldT = -oldT; }
            return new[] { oldR, oldS, oldT };
        }

        private static ulong MultiplyModulo(ulong a, ulong b, ulong modulus)
        {
            if (a <= uint.MaxValue && b <= uint.MaxValue) return a * b % modulus;
            return (ulong)((BigInteger)a * b % modulus);
        }

        private static ulong PowerModulo(ulong value, ulong exponent, ulong modulus)
        {
            ulong result = 1;
            while (exponent != 0)
            {
                if ((exponent & 1) != 0) result = MultiplyModulo(result, value, modulus);
                exponent >>= 1;
                if (exponent != 0) value = MultiplyModulo(value, value, modulus);
            }
            return result;
        }

        // Sorenson and Webster, https://arxiv.org/abs/1509.00864:
        // the first composite passing all twelve bases exceeds the entire Int64 range.
        private static readonly uint[] PrimalityBases = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };

        private static bool IsPrimeUnsigned(ulong value)
        {
            if (value < 2) return false;
            foreach (uint prime in PrimalityBases)
            {
                if (value % prime == 0) return value == prime;
            }
            ulong oddPart = value - 1;
            int powersOfTwo = 0;
            while ((oddPart & 1) == 0) { oddPart >>= 1; powersOfTwo++; }
            foreach (uint witness in PrimalityBases)
            {
                ulong residue = PowerModulo(witness, oddPart, value);
                if (residue == 1 || residue == value - 1) continue;
                bool passed = false;
                for (int i = 1; i < powersOfTwo; i++)
                {
                    residue = MultiplyModulo(residue, residue, value);
                    if (residue == value - 1) { passed = true; break; }
                }
                if (!passed) return false;
            }
            return true;
        }

        private static void FactorInteger(ulong value, List<ulong> factors)
        {
            if (value == 1) return;
            if (IsPrimeUnsigned(value)) { factors.Add(value); return; }
            ulong divisor = FindDivisor(value);
            FactorInteger(divisor, factors);
            FactorInteger(value / divisor, factors);
        }

        private static ulong FindDivisor(ulong value)
        {
            foreach (uint prime in PrimalityBases)
                if (value % prime == 0) return prime;

            // Brent's batched Pollard rho. Restart failed cycles; a failed split is not primality evidence.
            // Bounded attempts plus exact trial division guarantee a finite fallback.
            for (ulong constant = 1; constant <= 32; constant++)
            {
                ulong y = 2, x = 0, saved = 0, divisor = 1;
                for (int length = 1; length <= 131072 && divisor == 1; length *= 2)
                {
                    x = y;
                    for (int i = 0; i < length; i++) y = (MultiplyModulo(y, y, value) + constant) % value;
                    for (int start = 0; start < length && divisor == 1; start += 64)
                    {
                        saved = y;
                        ulong product = 1;
                        int count = Math.Min(64, length - start);
                        for (int i = 0; i < count; i++)
                        {
                            y = (MultiplyModulo(y, y, value) + constant) % value;
                            ulong difference = x > y ? x - y : y - x;
                            product = MultiplyModulo(product, difference, value);
                        }
                        divisor = UnsignedGcd(product, value);
                    }
                }
                if (divisor == value)
                {
                    for (int i = 0; i < 64; i++)
                    {
                        saved = (MultiplyModulo(saved, saved, value) + constant) % value;
                        divisor = UnsignedGcd(x > saved ? x - saved : saved - x, value);
                        if (divisor != 1) break;
                    }
                }
                if (divisor > 1 && divisor < value) return divisor;
            }
            for (ulong divisor = 41; divisor <= value / divisor; divisor += 2)
                if (value % divisor == 0) return divisor;
            return value;
        }

        private static long AccumulateDigits(int[] digits, int radix, bool leastSignificantFirst)
        {
            if (digits == null) throw new ArgumentNullException(nameof(digits));
            if (radix < 2) throw new ArgumentOutOfRangeException(nameof(radix));
            long result = 0;
            for (int i = 0; i < digits.Length; i++)
            {
                int digit = digits[leastSignificantFirst ? digits.Length - i - 1 : i];
                if (digit < 0 || digit >= radix) throw new ArgumentOutOfRangeException(nameof(digits), "Digits must be in the range [0, radix).");
                result = checked(result * radix + digit);
            }
            return result;
        }
    }
}
