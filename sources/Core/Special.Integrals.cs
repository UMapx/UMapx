using System;
using System.Numerics;

namespace UMapx.Core
{
    public static partial class Special
    {
        private static Complex ExponentialIntegralE1(Complex z)
        {
            if (z.Magnitude > 4) return Complex.Exp(-z) * GammaFraction(Complex.Zero, z);
            Complex term = -z, sum = term;
            for (int k = 2; k < 1000; k++)
            {
                term *= -z / k;
                Complex add = term / k;
                sum += add;
                if (add.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
            }
            return -EulerGamma - Complex.Log(z) - sum;
        }

        private static Complex ExponentialIntegral(Complex z)
        {
            if (z == Complex.Zero) return new Complex(double.NegativeInfinity, 0);
            if (!IsFinite(z)) return ComplexNaN;
            if (z.Real < 0)
                return -ExponentialIntegralE1(-z) + new Complex(0, z.Imaginary > 0 ? Math.PI : z.Imaginary < 0 ? -Math.PI : 0);
            if (z.Magnitude >= 20)
            {
                Complex term = 1, sum = 1;
                double previous = double.PositiveInfinity;
                for (int k = 1; k < 1000; k++)
                {
                    term *= k / z;
                    if (term.Magnitude >= previous) break;
                    sum += term;
                    if (term.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
                    previous = term.Magnitude;
                }
                return Complex.Exp(z) / z * sum + new Complex(0, z.Imaginary > 0 ? Math.PI : z.Imaginary < 0 ? -Math.PI : 0);
            }
            Complex power = z, series = power;
            for (int k = 2; k < 1000; k++)
            {
                power *= z / k;
                Complex add = power / k;
                series += add;
                if (add.Magnitude <= SpecialEpsilon * series.Magnitude) break;
            }
            return EulerGamma + Complex.Log(z) + series;
        }

        private static Complex TrigonometricIntegral(Complex z, bool sine)
        {
            if (z == Complex.Zero) return sine ? Complex.Zero : new Complex(double.NegativeInfinity, 0);
            if (!IsFinite(z)) return ComplexNaN;
            if (z.Real < 0) return sine ? -TrigonometricIntegral(-z, true)
                : TrigonometricIntegral(-z, false) + new Complex(0, z.Imaginary < 0 ? -Math.PI : Math.PI);
            if (z.Magnitude < 1)
            {
                Complex term = sine ? z : -z * z / 4, sum = term;
                for (int k = 1; k < 100; k++)
                {
                    term *= sine ? -z * z * (2 * k - 1) / ((2.0 * k + 1) * (2 * k + 1) * (2 * k))
                        : -z * z * (2 * k) / ((2.0 * k + 2) * (2 * k + 2) * (2 * k + 1));
                    sum += term;
                    if (term.Magnitude <= SpecialEpsilon * sum.Magnitude) break;
                }
                return sine ? sum : EulerGamma + Complex.Log(z) + sum;
            }
            Complex positive = ExponentialIntegral(Complex.ImaginaryOne * z);
            Complex negative = ExponentialIntegral(-Complex.ImaginaryOne * z);
            if (z.Real == 0)
            {
                if (z.Imaginary > 0) positive += Complex.ImaginaryOne * Math.PI;
                else negative -= Complex.ImaginaryOne * Math.PI;
            }
            return sine ? (positive - negative) / (2 * Complex.ImaginaryOne) - Math.PI / 2 : (positive + negative) / 2;
        }

        private static Complex OwenValue(Complex h, Complex a)
        {
            if (!IsFinite(h) || !IsFinite(a)) return ComplexNaN;
            if (a == Complex.Zero) return Complex.Zero;
            Complex angle = Complex.Atan(a);
            if (h == Complex.Zero) return angle / (2 * Math.PI);
            Complex previous = ComplexNaN;
            for (int panels = 2; panels <= 1024; panels *= 2)
            {
                Complex value = angle / (2 * Math.PI) * IntegrateGauss(t =>
                {
                    Complex cosine = Complex.Cos(angle * t);
                    return Complex.Exp(-h * h / (2 * cosine * cosine));
                }, 1, panels);
                if ((value - previous).Magnitude <= 1e-13 * value.Magnitude) return value;
                previous = value;
            }
            return previous;
        }
    }
}
