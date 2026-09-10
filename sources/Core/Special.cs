using System;
using System.Numerics;

namespace UMapx.Core
{
    /// <summary>
    /// Used to implement special mathematical functions.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Special_functions
    /// </remarks>
    public static partial class Special
    {
        #region Chebyshev polynomial
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static float ChebyshevT(float x, int n)
        {
            return (float)ChebyshevValue((Complex)x, n, false).Real;
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static Complex32 ChebyshevT(Complex32 x, int n)
        {
            return (Complex32)ChebyshevValue((Complex)x, n, false);
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static float ChebyshevU(float x, int n)
        {
            return (float)ChebyshevValue((Complex)x, n, true).Real;
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static Complex32 ChebyshevU(Complex32 x, int n)
        {
            return (Complex32)ChebyshevValue((Complex)x, n, true);
        }
        #endregion

        #region Abel polynomial
        /// <summary>
        /// Returns the value of the Abel polynomial.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Power</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static float Abel(float x, float a, int n)
        {
            if (n < 0) return float.NaN;
            if (n == 0) return 1;
            return (float)((Complex)x * Complex.Pow((Complex)x - n * (Complex)a, n - 1)).Real;
        }
        /// <summary>
        /// Returns the value of the Abel polynomial.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Complex power</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Abel(Complex32 x, Complex32 a, int n)
        {
            if (n < 0) return Complex32.NaN;
            if (n == 0) return 1;
            return (Complex32)((Complex)x * Complex.Pow((Complex)x - n * (Complex)a, n - 1));
        }
        #endregion

        #region Laguerre polynomial
        /// <summary>
        /// Returns the value of the Laguerre polynomial.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Power</param>
        /// <param name="k">Order</param>
        /// <returns>Value</returns>
        public static float Laguerre(float x, float a, int k)
        {
            return (float)OrthogonalPolynomial((Complex)x, (Complex)a, k, 0).Real;
        }
        /// <summary>
        /// Returns the value of the Laguerre polynomial.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Power</param>
        /// <param name="k">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Laguerre(Complex32 x, Complex32 a, int k)
        {
            return (Complex32)OrthogonalPolynomial((Complex)x, (Complex)a, k, 0);
        }
        #endregion

        #region Legendre polynomial
        /// <summary>
        /// Returns the value of the Legendre polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="m">Order</param>
        /// <returns>Value</returns>
        public static float Legendre(float x, int m)
        {
            return (float)OrthogonalPolynomial((Complex)x, Complex.Zero, m, 1).Real;
        }
        /// <summary>
        /// Returns the value of the Legendre polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="m">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Legendre(Complex32 x, int m)
        {
            return (Complex32)OrthogonalPolynomial((Complex)x, Complex.Zero, m, 1);
        }
        #endregion

        #region Hermite polynomial
        /// <summary>
        /// Returns the value of the Hermite polynomial.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="m">Order</param>
        /// <returns>Value</returns>
        public static float Hermite(float x, int m)
        {
            return (float)OrthogonalPolynomial((Complex)x, Complex.Zero, m, 2).Real;
        }
        /// <summary>
        /// Returns the value of the Hermite polynomial.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="m">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Hermite(Complex32 x, int m)
        {
            return (Complex32)OrthogonalPolynomial((Complex)x, Complex.Zero, m, 2);
        }
        #endregion

        #region Gegenbauer polynomial
        /// <summary>
        /// Returns the value of the Gegenbauer polynomial.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Power</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static float Gegenbauer(float x, float a, int n)
        {
            return (float)OrthogonalPolynomial((Complex)x, (Complex)a, n, 3).Real;
        }
        /// <summary>
        /// Returns the value of the Gegenbauer polynomial.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Power</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Gegenbauer(Complex32 x, Complex32 a, int n)
        {
            return (Complex32)OrthogonalPolynomial((Complex)x, (Complex)a, n, 3);
        }
        #endregion

        #region Sinc function
        /// <summary>
        /// Returns the value of the normalized cardinal sine function: f(x) = sin(πx) / (πx).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Sinc(float x)
        {
            return Special.Sinc(x, Maths.Pi);
        }
        /// <summary>
        /// Returns the value of the normalized cardinal sine function: f(x) = sin(πx) / (πx).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Sinc(Complex32 x)
        {
            return Special.Sinc(x, Maths.Pi);
        }
        /// <summary>
        /// Returns the value of the cardinal sine function with the parameter: f(x, a) = sin(ax) / (ax).
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Sinc(float x, float a)
        {
            var ax = a * x;

            if (ax == 0)
                return 1;

            return Maths.Sin(ax) / ax;
        }
        /// <summary>
        /// Returns the value of the cardinal sine function with the parameter: f(x, a) = sin(ax) / (ax).
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Sinc(Complex32 x, Complex32 a)
        {
            var ax = a * x;

            if (ax == Complex32.Zero)
                return Complex32.One;

            return Maths.Sin(ax) / ax;
        }
        #endregion

        #region Guderman & Hartley functions
        /// <summary>
        /// Returns the value of the inverse Guderman function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Agd(float x)
        {
            // gd^{-1}(x) = artanh(sin(x))
            return Maths.Atanh(Maths.Sin(x));
        }
        /// <summary>
        /// Returns the value of the inverse Guderman function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Agd(Complex32 x)
        {
            // gd^{-1}(x) = artanh(sin(x))
            return Maths.Atanh(Maths.Sin(x));
        }
        /// <summary>
        /// Returns the value of the Guderman function.
        /// </summary>
        /// <param name="x">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Gd(float x)
        {
            return (float)(2 * Math.Atan(Math.Tanh((double)x / 2)));
        }
        /// <summary>
        /// Returns the value of the Guderman function.
        /// </summary>
        /// <param name="x">Angle in radians</param>
        /// <returns>Value</returns>
        public static Complex32 Gd(Complex32 x)
        {
            return (Complex32)(2 * Complex.Atan(2 * LogisticValue((Complex)x) - 1));
        }
        /// <summary>
        /// Returns the value of the function Cas(x).
        /// </summary>
        /// <param name="theta">Theta</param>
        /// <returns>Value</returns>
        public static float Cas(float theta)
        {
            return Maths.Cos(theta) + Maths.Sin(theta);
        }
        /// <summary>
        /// Returns the value of the function Cas(x).
        /// </summary>
        /// <param name="theta">Theta</param>
        /// <returns>Value</returns>
        public static Complex32 Cas(Complex32 theta)
        {
            return Maths.Cos(theta) + Maths.Sin(theta);
        }
        #endregion

        #region Rademacher function
        /// <summary>
        /// Returns the value of the Radamecher function.
        /// </summary>
        /// <param name="t">Value [0, 1]</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static float Rademacher(float t, int n)
        {
            return (float)RademacherValue((Complex)t, n).Real;
        }
        /// <summary>
        /// Returns the value of the Radamecher function.
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Rademacher(Complex32 z, int n)
        {
            return (Complex32)RademacherValue((Complex)z, n);
        }
        #endregion

        #region Heavyside delta-function
        /// <summary>
        /// Returns the value of the Heaviside delta function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="k">Smoothing factor</param>
        /// <returns>Value</returns>
        public static float Heaviside(float x, float k)
        {
            return (float)LogisticValue(2 * (Complex)k * (Complex)x).Real;
        }
        /// <summary>
        /// Returns the value of the Heaviside delta function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="k">Smoothing factor</param>
        /// <returns>Value</returns>
        public static Complex32 Heaviside(Complex32 x, Complex32 k)
        {
            return (Complex32)LogisticValue(2 * (Complex)k * (Complex)x);
        }
        #endregion

        #region Mahler function
        /// <summary>
        /// Returns the value of the Mahler function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="t">Value</param>
        /// <returns>Value</returns>
        public static float Mahler(float x, float t)
        {
            return Maths.Exp(x * (1.0f + t - Maths.Pow(Maths.E, t)));
        }
        /// <summary>
        /// Returns the value of the Mahler function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="t">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Mahler(Complex32 x, Complex32 t)
        {
            return Maths.Exp(x * (1.0f + t - Maths.Pow(Maths.E, t)));
        }
        #endregion

        #region Gompertz function
        /// <summary>
        /// Gets the value of the Gompertz function.
        /// </summary>
        /// <param name="t">Value</param>
        /// <param name="a">Upper asymptote</param>
        /// <param name="b">Growth parameter</param>
        /// <param name="c">Growth rate</param>
        /// <returns>Value</returns>
        public static float Gompertz(float t, float a, float b, float c)
        {
            return a * Maths.Exp(-b * Maths.Exp(-c * t));
        }
        /// <summary>
        /// Gets the value of the Gompertz function.
        /// </summary>
        /// <param name="t">Value</param>
        /// <param name="a">Upper asymptote</param>
        /// <param name="b">Growth parameter</param>
        /// <param name="c">Growth rate</param>
        /// <returns>Value</returns>
        public static Complex32 Gompertz(Complex32 t, Complex32 a, Complex32 b, Complex32 c)
        {
            return a * Maths.Exp(-b * Maths.Exp(-c * t));
        }
        #endregion

        #region Dirac delta-function
        /// <summary>
        /// Returns the value of the Dirac delta function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Coefficient</param>
        /// <returns>Value</returns>
        public static float Dirac(float x, float a)
        {
            float s = Maths.Sqrt(Maths.Pi);
            float b = 1.0f / Math.Abs(a) / s;
            float c = Maths.Pow(x / a, 2);
            float e = Maths.Exp(-c);
            return b * e;
        }
        /// <summary>
        /// Returns the value of the Dirac delta function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Coefficient</param>
        /// <returns>Value</returns>
        public static Complex32 Dirac(Complex32 x, Complex32 a)
        {
            Complex32 s = Maths.Sqrt(Maths.Pi);
            Complex32 b = 1.0f / Maths.Abs(a) / s;
            Complex32 c = Maths.Pow(x / a, 2);
            Complex32 e = Maths.Exp(-c);
            return b * e;
        }
        #endregion

        #region Logistic function
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Lower asymptote</param>
        /// <param name="k">Upper asymptote</param>
        /// <param name="b">Growth rate</param>
        /// <param name="v">Affect</param>
        /// <param name="q">Central moment</param>
        /// <param name="c">Offset</param>
        /// <returns>Value</returns>
        public static float Logistic(float x, float a, float k, float b, float v, float q, float c)
        {
            return a + (k - a) / Maths.Pow(c + q * Maths.Exp(-b * x), 1.0f / v);
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Lower asymptote</param>
        /// <param name="k">Upper asymptote</param>
        /// <param name="b">Growth rate</param>
        /// <param name="v">Affect</param>
        /// <param name="q">Central moment</param>
        /// <param name="c">Offset</param>
        /// <returns>Value</returns>
        public static Complex32 Logistic(Complex32 x, Complex32 a, Complex32 k, Complex32 b, Complex32 v, Complex32 q, Complex32 c)
        {
            return a + (k - a) / Maths.Pow(c + q * Maths.Exp(-b * x), 1.0f / v);
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Lower asymptote</param>
        /// <param name="k">Upper asymptote</param>
        /// <param name="b">Growth rate</param>
        /// <returns>Value</returns>
        public static float Logistic(float x, float a, float k, float b)
        {
            return (float)((Complex)a + ((Complex)k - (Complex)a) * LogisticValue((Complex)b * (Complex)x)).Real;
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Lower asymptote</param>
        /// <param name="k">Upper asymptote</param>
        /// <param name="b">Growth rate</param>
        /// <returns>Value</returns>
        public static Complex32 Logistic(Complex32 x, Complex32 a, Complex32 k, Complex32 b)
        {
            return (Complex32)((Complex)a + ((Complex)k - (Complex)a) * LogisticValue((Complex)b * (Complex)x));
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Logistic(float x)
        {
            return (float)LogisticValue((Complex)x).Real;
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Logistic(Complex32 x)
        {
            return (Complex32)LogisticValue((Complex)x);
        }
        #endregion

        #region Elrang B and C functions
        /// <summary>
        /// Returns the value of the Erlang C-function.
        /// </summary>
        /// <param name="y">First parameter</param>
        /// <param name="v">Second parameter</param>
        /// <param name="t">Time parameter</param>
        /// <returns>Value</returns>
        public static float Erlang(float y, int v, float t)
        {
            Complex traffic = (Complex)y, blocking = ErlangBlocking(traffic, v);
            return (float)(v * blocking / (v - traffic + traffic * blocking) * Complex.Exp(-(v - traffic) * (Complex)t)).Real;
        }
        /// <summary>
        /// Returns the value of the Erlang C-function.
        /// </summary>
        /// <param name="y">First parameter</param>
        /// <param name="v">Second parameter</param>
        /// <param name="t">Time parameter</param>
        /// <returns>Value</returns>
        public static Complex32 Erlang(Complex32 y, int v, Complex32 t)
        {
            Complex traffic = (Complex)y, blocking = ErlangBlocking(traffic, v);
            return (Complex32)(v * blocking / (v - traffic + traffic * blocking) * Complex.Exp(-(v - traffic) * (Complex)t));
        }
        /// <summary>
        /// Returns the value of the Erlang B-function.
        /// </summary>
        /// <param name="y">First parameter</param>
        /// <param name="v">Second parameter</param>
        /// <returns>Value</returns>
        public static float Erlang(float y, int v)
        {
            return (float)ErlangBlocking((Complex)y, v).Real;
        }
        /// <summary>
        /// Returns the value of the Erlang B-function.
        /// </summary>
        /// <param name="y">First parameter</param>
        /// <param name="v">Second parameter</param>
        /// <returns>Value</returns>
        public static Complex32 Erlang(Complex32 y, int v)
        {
            return (Complex32)ErlangBlocking((Complex)y, v);
        }
        #endregion

        #region Lambert W-function
        /// <summary>
        /// Returns the value of the Lambert W-function.
        /// </summary>
        /// <param name="x">Value [-1/e,+inf)</param>
        /// <param name="k">Branch</param>
        /// <returns>Value</returns>
        public static float LambertW(float x, int k = 0)
        {
            if ((k != 0 && k != -1) || x < -0.36787944117144233f || (k == -1 && x >= 0)) return float.NaN;
            if (x == -0.36787944117144233f) return -1;
            if (float.IsPositiveInfinity(x)) return float.PositiveInfinity;
            return (float)LambertValue((Complex)(double)x, k).Real;
        }
        /// <summary>
        /// Returns the value of the Lambert W-function.
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="k">Branch</param>
        /// <returns>Value</returns>
        public static Complex32 LambertW(Complex32 z, int k = 0)
        {
            return (Complex32)LambertValue((Complex)z, k);
        }
        /// <summary>
        /// Returns the value of the square super-root.
        /// </summary>
        /// <param name="x">Value [1,+inf)</param>
        /// <param name="k">Branch</param>
        /// <returns>Value</returns>
        public static float Ssqrt(float x, int k = 0)
        {
            if (x <= 0 || (k != 0 && k != -1)) return float.NaN;
            if (x == 1 && k == 0) return 1;
            Complex value = LambertValue((Complex)Math.Log(x), k);
            return Math.Abs(value.Imaginary) <= 1e-12 ? (float)Math.Exp(value.Real) : float.NaN;
        }
        /// <summary>
        /// Returns the value of the square super-root.
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="k">Branch</param>
        /// <returns>Value</returns>
        public static Complex32 Ssqrt(Complex32 z, int k = 0)
        {
            if (z == 1 && k == 0) return 1;
            return (Complex32)Complex.Exp(LambertValue(Complex.Log((Complex)z), k));
        }
        #endregion

        #region Fresnel integral functions
        /// <summary>
        /// Returns the value of the Fresnel integral C(x).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Fresnelc(float x)
        {
            return (float)FresnelValue((Complex)x, false).Real;
        }
        /// <summary>
        /// Returns the value of the Fresnel integral C(x).
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Fresnelc(Complex32 z)
        {
            return (Complex32)FresnelValue((Complex)z, false);
        }
        /// <summary>
        /// Returns the value of the Fresnel integral S(x).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Fresnels(float x)
        {
            return (float)FresnelValue((Complex)x, true).Real;
        }
        /// <summary>
        /// Returns the value of the Fresnel integral S(x).
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Fresnels(Complex32 z)
        {
            return (Complex32)FresnelValue((Complex)z, true);
        }
        #endregion

        #region Owen's T-function
        /// <summary>
        /// Returns the value of the Owen T function.
        /// </summary>
        /// <param name="h">First value</param>
        /// <param name="a">Second value</param>
        /// <returns>Value</returns>
        public static float Owen(float h, float a)
        {
            return (float)OwenValue((Complex)h, (Complex)a).Real;
        }
        /// <summary>
        /// Returns the value of the Owen T function.
        /// </summary>
        /// <param name="h">First value</param>
        /// <param name="a">Second value</param>
        /// <returns>Value</returns>
        public static Complex32 Owen(Complex32 h, Complex32 a)
        {
            return (Complex32)OwenValue((Complex)h, (Complex)a);
        }

        #endregion

        #region Riemann's Zeta function
        /// <summary>
        /// Returns the value of the Riemann zeta ζ(s) on the principal branch (real s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <returns>Value</returns>
        public static float Zeta(float s)
        {
            if (s == 1) return float.PositiveInfinity;
            if (float.IsPositiveInfinity(s)) return 1;
            return (float)ZetaValue((Complex)s).Real;
        }
        /// <summary>
        /// Returns the value of the Riemann zeta ζ(s) on the principal branch (complex s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Zeta(Complex32 s)
        {
            return (Complex32)ZetaValue((Complex)s);
        }



        #endregion



        #region Gamma functions

        /// <summary>
        /// Returns the value of the Euler Gamma function: Gamma(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Gamma(float x)
        {
            return (float)GammaValue((double)x);
        }
        /// <summary>
        /// Returns the value of the Euler Gamma function: Gamma(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Gamma(Complex32 x)
        {
            return (Complex32)GammaValue((Complex)x);
        }

        /// <summary>
        /// Returns log(abs(Gamma(x))) for real arguments. Gamma poles return NaN.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float LogGamma(float x)
        {
            return (float)GammaLog((double)x);
        }
        /// <summary>
        /// Returns analytic log-gamma with its cut on the negative real axis.
        /// Its imaginary part is not reduced modulo 2*pi; the upper side is used on the cut.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 LogGamma(Complex32 z)
        {
            return (Complex32)GammaLog((Complex)z);
        }

        /// <summary>
        /// Returns the value of the Digamma function: ψ(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float DiGamma(float x)
        {
            if (float.IsPositiveInfinity(x)) return float.PositiveInfinity;
            return (float)Polygamma((Complex)x, false).Real;
        }
        /// <summary>
        /// Returns the value of the Digamma function: ψ(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 DiGamma(Complex32 x)
        {
            return (Complex32)Polygamma((Complex)x, false);
        }

        /// <summary>
        /// Returns the value of the Trigamma function: ψ1(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float TriGamma(float x)
        {
            if (float.IsPositiveInfinity(x)) return 0;
            return (float)Polygamma((Complex)x, true).Real;
        }
        /// <summary>
        /// Returns the value of the Trigamma function: ψ1(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 TriGamma(Complex32 x)
        {
            return (Complex32)Polygamma((Complex)x, true);
        }

        /// <summary>
        /// Returns the value of the incomplete upper Gamma function: Q(s, x) = Γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float GammaQ(float s, float x)
        {
            return (float)IncompleteGamma((double)s, (double)x, true, true);
        }
        /// <summary>
        /// Returns the value of the incomplete upper Gamma function: Q(s, x) = Γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 GammaQ(Complex32 s, Complex32 x)
        {
            return (Complex32)IncompleteGamma((Complex)s, (Complex)x, true, true);
        }

        /// <summary>
        /// Returns the value of an incomplete lower Gamma function: P(s, x) = γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float GammaP(float s, float x)
        {
            return (float)IncompleteGamma((double)s, (double)x, false, true);
        }
        /// <summary>
        /// Returns the value of an incomplete lower Gamma function: P(s, x) = γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 GammaP(Complex32 s, Complex32 x)
        {
            return (Complex32)IncompleteGamma((Complex)s, (Complex)x, false, true);
        }

        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float GammaIncomplete(float s, float x)
        {
            return (float)IncompleteGamma((double)s, (double)x, false, false);
        }
        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 GammaIncomplete(Complex32 s, Complex32 x)
        {
            return (Complex32)IncompleteGamma((Complex)s, (Complex)x, false, false);
        }

        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x) (complemented).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float GammaIncompleteComplemented(float s, float x)
        {
            return (float)IncompleteGamma((double)s, (double)x, true, false);
        }
        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x) (complemented).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 GammaIncompleteComplemented(Complex32 s, Complex32 x)
        {
            return (Complex32)IncompleteGamma((Complex)s, (Complex)x, true, false);
        }



        #endregion

        #region Generalized error function
        /// <summary>
        /// Returns n!/sqrt(pi) times the integral of exp(-t^n) from zero to x.
        /// </summary>
        /// <param name="x">Real integration endpoint</param>
        /// <param name="n">Order [0, +inf)</param>
        /// <returns>Value</returns>
        public static float Gerf(float x, int n)
        {
            return (float)GeneralizedErf((Complex)x, n).Real;
        }
        /// <summary>
        /// Returns the entire continuation of n!/sqrt(pi) times the integral of exp(-t^n).
        /// </summary>
        /// <param name="x">Complex integration endpoint</param>
        /// <param name="n">Order [0, +inf)</param>
        /// <returns>Value</returns>
        public static Complex32 Gerf(Complex32 x, int n)
        {
            return (Complex32)GeneralizedErf((Complex)x, n);
        }
        /// <summary>
        /// Returns the value of the generalized error function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Gerf(float x)
        {
            return Gerf(x, 2);
        }
        /// <summary>
        /// Returns the value of the generalized error function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Gerf(Complex32 x) => Gerf(x, 2);
        #endregion

        #region Factorial function
        /// <summary>
        /// Returns the natural logarithm of the factorial of a number log(n!).
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Value</returns>
        public static float LogFactorial(float n)
        {
            return (float)GammaLog((double)n + 1);
        }
        /// <summary>
        /// Returns the natural logarithm of the factorial of a number log(n!).
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 LogFactorial(Complex32 z)
        {
            return (Complex32)GammaLog((Complex)z + 1);
        }
        /// <summary>
        /// Returns the factorial of a number.
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Value</returns>
        public static double Factorial(float n)
        {
            if (n >= 0 && n <= 170 && n == Math.Floor(n)) return A000142[(int)n];
            return GammaValue((double)n + 1);
        }
        /// <summary>
        /// Returns the factorial of a number.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Factorial(Complex32 z)
        {
            return (Complex32)GammaValue((Complex)z + 1);
        }
        /// <summary>
        /// Returns the decreasing factorial of a number.
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static double FactorialDown(float n, float k)
        {
            return FactorialProduct((Complex)n, (Complex)k, false).Real;
        }
        /// <summary>
        /// Returns the decreasing factorial of a number.
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static Complex32 FactorialDown(Complex32 z, Complex32 k)
        {
            return (Complex32)FactorialProduct((Complex)z, (Complex)k, false);
        }
        /// <summary>
        /// Returns the increasing factorial of a number (Pochhammer symbol).
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static float FactorialUp(float n, float k)
        {
            return (float)FactorialProduct((Complex)n, (Complex)k, true).Real;
        }
        /// <summary>
        /// Returns the increasing factorial of a number (Pochhammer symbol).
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static Complex32 FactorialUp(Complex32 z, Complex32 k)
        {
            return (Complex32)FactorialProduct((Complex)z, (Complex)k, true);
        }
        #endregion

        #region Binomial function
        /// <summary>
        /// Returns the value of binomial coefficients: C(n, k) = n! / k! / (n-k)! for k > 0.
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static double Binomial(float n, float k)
        {
            return BinomialValue((Complex)n, (Complex)k).Real;
        }
        /// <summary>
        /// Returns the value of binomial coefficients: C(n, k) = n! / k! / (n-k)! for k > 0.
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Binomial(Complex32 n, Complex32 k)
        {
            return (Complex32)BinomialValue((Complex)n, (Complex)k);
        }
        /// <summary>
        /// Returns the natural logarithm of binomial coefficients: log(C(n, k)) = log(n!) - log(k!) - log(n-k!).
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static float LogBinomial(float n, float k)
        {
            if (k < 0 || (n >= 0 && n == Math.Floor(n) && k > n)) return float.NegativeInfinity;
            return (float)(GammaLog((double)n + 1) - GammaLog((double)k + 1) - GammaLog((double)n - k + 1));
        }
        /// <summary>
        /// Returns the natural logarithm of binomial coefficients: log(C(n, k)) = log(n!) - log(k!) - log(n-k!).
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static Complex32 LogBinomial(Complex32 n, Complex32 k)
        {
            return (Complex32)(GammaLog((Complex)n + 1) - GammaLog((Complex)k + 1) - GammaLog((Complex)n - (Complex)k + 1));
        }
        #endregion

        #region Laplace functions
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="inverse">Reverse function or not</param>
        /// <returns>Value</returns>
        public static float Erf(float x, bool inverse)
        {
            return (float)(inverse ? InverseErf((double)x) : ErfValue((double)x));
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="inverse">Reverse function or not</param>
        /// <returns>Value</returns>
        public static Complex32 Erf(Complex32 x, bool inverse)
        {
            return (Complex32)(inverse ? InverseErf((Complex)x) : ErfValue((Complex)x));
        }
        /// <summary>
        /// Returns the value of the imaginary error function.
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static float Erfi(float x)
        {
            if (Math.Abs(x) > 27) return x < 0 ? float.NegativeInfinity : float.PositiveInfinity;
            return (float)(-Complex.ImaginaryOne * ErfValue(Complex.ImaginaryOne * (double)x)).Real;
        }
        /// <summary>
        /// Returns the value of the imaginary error function.
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static Complex32 Erfi(Complex32 x)
        {
            return (Complex32)(-Complex.ImaginaryOne * ErfValue(Complex.ImaginaryOne * (Complex)x));
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static float Erf(float x)
        {
            return Erf(x, false);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static Complex32 Erf(Complex32 x)
        {
            return Erf(x, false);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="a">The lower boundary of the normalization</param>
        /// <param name="b">The upper limit of the normalization</param>
        /// <returns>Value</returns>
        public static float Erf(float x, float a, float b)
        {
            return (float)ErfValue(((double)x - (double)a) / (double)b);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="a">The lower boundary of the normalization</param>
        /// <param name="b">The upper limit of the normalization</param>
        /// <returns>Value</returns>
        public static Complex32 Erf(Complex32 x, Complex32 a, Complex32 b)
        {
            return (Complex32)ErfValue(((Complex)x - (Complex)a) / (Complex)b);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static float Erfc(float x)
        {
            return (float)ErfcValue((double)x);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static Complex32 Erfc(Complex32 x)
        {
            return (Complex32)ErfcValue((Complex)x);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="a">The lower boundary of the normalization</param>
        /// <param name="b">The upper limit of the normalization</param>
        /// <returns>Value</returns>
        public static float Erfc(float x, float a, float b)
        {
            return (float)ErfcValue(((double)x - (double)a) / (double)b);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="a">The lower boundary of the normalization</param>
        /// <param name="b">The upper limit of the normalization</param>
        /// <returns>Value</returns>
        public static Complex32 Erfc(Complex32 x, Complex32 a, Complex32 b)
        {
            return (Complex32)ErfcValue(((Complex)x - (Complex)a) / (Complex)b);
        }



        #endregion

        #region Dawson function
        /// <summary>
        /// Returns the value of the D- / D + Dawson function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="positive">D- or D+</param>
        /// <returns>Value</returns>
        public static float Dawson(float x, bool positive)
        {
            if (!positive && Math.Abs(x) > 27) return x < 0 ? float.NegativeInfinity : float.PositiveInfinity;
            return (float)(positive ? DawsonValue((Complex)x) : -Complex.ImaginaryOne * DawsonValue(Complex.ImaginaryOne * (double)x)).Real;
        }
        /// <summary>
        /// Returns the value of the D- / D + Dawson function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="positive">D- or D+</param>
        /// <returns>Value</returns>
        public static Complex32 Dawson(Complex32 x, bool positive)
        {
            return (Complex32)(positive ? DawsonValue((Complex)x) : -Complex.ImaginaryOne * DawsonValue(Complex.ImaginaryOne * (Complex)x));
        }
        #endregion

        #region Faddeeva function
        /// <summary>
        /// Returns the value of the Faddeeva function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Faddeeva(float x)
        {
            return (Complex32)FaddeevaValue((Complex)x);
        }
        /// <summary>
        /// Returns the value of the Faddeeva function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Faddeeva(Complex32 x)
        {
            return (Complex32)FaddeevaValue((Complex)x);
        }
        #endregion

        #region Q-function
        /// <summary>
        /// Returns the value of a Q function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="inverse">Inverse function or not</param>
        /// <returns>Value</returns>
        public static float Q(float x, bool inverse = false)
        {
            return (float)(inverse ? Math.Sqrt(2) * InverseErfc(2.0 * x) : 0.5 * ErfcValue((double)x / Math.Sqrt(2)));
        }
        /// <summary>
        /// Returns the value of a Q function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="inverse">Inverse function or not</param>
        /// <returns>Value</returns>
        public static Complex32 Q(Complex32 x, bool inverse = false)
        {
            if (inverse && x.Imag == 0 && x.Real >= 0 && x.Real <= 1) return new Complex32((float)(Math.Sqrt(2) * InverseErfc(2.0 * x.Real)), 0);
            return (Complex32)(inverse ? Math.Sqrt(2) * InverseErf(1 - 2 * (Complex)x) : 0.5 * ErfcValue((Complex)x / Math.Sqrt(2)));
        }
        #endregion

        #region Hypergeometric function
        /// <summary>
        /// Returns the value of a hypergeometric function.
        /// </summary>
        /// <remarks>
        /// This version of the hypergeometric function is found in the Russian literature and is indicated: F(a,b,c,z).
        /// More information can be found on the website:
        /// https://en.wikipedia.org/wiki/Hypergeometric_function
        /// </remarks>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <param name="c">Value</param>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static float Hypergeom(float a, float b, float c, float z)
        {
            Complex value = Hypergeometric2F1((Complex)a, (Complex)b, (Complex)c, (Complex)z);
            return Math.Abs(value.Imaginary) <= 1e-12 * (1 + Math.Abs(value.Real)) ? (float)value.Real : float.NaN;
        }
        /// <summary>
        /// Returns the value of a hypergeometric function.
        /// </summary>
        /// <remarks>
        /// This version of the hypergeometric function is found in the Russian literature and is indicated: F(a,b,c,z).
        /// More information can be found on the website:
        /// https://en.wikipedia.org/wiki/Hypergeometric_function
        /// </remarks>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <param name="c">Value</param>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Hypergeom(Complex32 a, Complex32 b, Complex32 c, Complex32 z)
        {
            return (Complex32)Hypergeometric2F1((Complex)a, (Complex)b, (Complex)c, (Complex)z);
        }
        /// <summary>
        /// Returns the value of a hypergeometric function.
        /// </summary>
        /// <remarks>
        /// The hypergeometric function can be used in several variations:
        /// F(a,b,z); F(a,~,z); F(~,b,z); F(~,~,z).
        /// Instead of the “~” sign, use the float.NaN value.
        /// More information can be found on the website:
        /// https://www.mathworks.com/help/symbolic/hypergeom.html#bt1nkmw-2
        /// </remarks>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static float Hypergeom(float a, float b, float z)
        {
            return (float)Hypergeometric1F1((Complex)a, (Complex)b, (Complex)z).Real;
        }
        /// <summary>
        /// Returns the value of a hypergeometric function.
        /// </summary>
        /// <remarks>
        /// The hypergeometric function can be used in several variations:
        /// F(a,b,z); F(a,~,z); F(~,b,z); F(~,~,z).
        /// Instead of the “~” sign, use the float.NaN value.
        /// More information can be found on the website:
        /// https://www.mathworks.com/help/symbolic/hypergeom.html#bt1nkmw-2
        /// </remarks>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Hypergeom(Complex32 a, Complex32 b, Complex32 z)
        {
            return (Complex32)Hypergeometric1F1((Complex)a, (Complex)b, (Complex)z);
        }
        #endregion

        #region Beta functions
        /// <summary>
        /// Returns the value of the beta function: B(a, b) = Gamma(a) * Gamma(b) / Gamma(a + b).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static float Beta(float a, float b)
        {
            return (float)BetaValue((double)a, (double)b);
        }
        /// <summary>
        /// Returns the value of the beta function: B(a, b) = Gamma(a) * Gamma(b) / Gamma(a + b).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Beta(Complex32 a, Complex32 b)
        {
            return (Complex32)Complex.Exp(BetaLog((Complex)a, (Complex)b));
        }
        /// <summary>
        /// Returns the value of the beta function: B(m, n) = (m - 1)! * (n - 1)! / (m + n - 1)!.
        /// </summary>
        /// <param name="m">Integer number</param>
        /// <param name="n">Integer number</param>
        /// <returns>Value</returns>
        public static double Beta(int m, int n)
        {
            return BetaValue((double)m, (double)n);
        }
        /// <summary>
        /// Returns the value of a derivative beta function: B'(a, b).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static float BetaDerivative(float a, float b)
        {
            return (float)(Complex.Exp(BetaLog((Complex)a, (Complex)b)) * (Polygamma((Complex)a, false) - Polygamma((Complex)a + (Complex)b, false))).Real;
        }
        /// <summary>
        /// Returns the value of a derivative beta function: B'(a, b).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static Complex32 BetaDerivative(Complex32 a, Complex32 b)
        {
            return (Complex32)(Complex.Exp(BetaLog((Complex)a, (Complex)b)) * (Polygamma((Complex)a, false) - Polygamma((Complex)a + (Complex)b, false)));
        }
        /// <summary>
        /// Returns the value of an incomplete beta function: Bx(a, b).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float BetaIncomplete(float a, float b, float x)
        {
            return (float)IncompleteBeta(a, b, x, false);
        }
        /// <summary>
        /// Returns the value of an incomplete beta function: Bx(a, b).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 BetaIncomplete(Complex32 a, Complex32 b, Complex32 x)
        {
            Complex aa = (Complex)a, bb = (Complex)b, z = (Complex)x;
            if (z == Complex.Zero) return Complex32.Zero;
            if (z == Complex.One) return (Complex32)Complex.Exp(BetaLog(aa, bb));
            return (Complex32)(Complex.Pow(z, aa) / aa * Hypergeometric2F1(aa, 1 - bb, aa + 1, z));
        }
        /// <summary>
        /// Returns the value of a regularized incomplete beta function: Ix(a, b).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float BetaIncompleteRegularized(float a, float b, float x)
        {
            return (float)IncompleteBeta(a, b, x, true);
        }
        /// <summary>
        /// Returns the value of a log-beta function.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static float LogBeta(float a, float b)
        {
            return (float)BetaLog((double)a, (double)b);
        }
        /// <summary>
        /// Returns the value of a log-beta function.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static Complex32 LogBeta(Complex32 a, Complex32 b)
        {
            return (Complex32)BetaLog((Complex)a, (Complex)b);
        }
        #endregion

        #region Integral functions
        /// <summary>
        /// Returns the value of the integral cosine.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Ci(float x)
        {
            if (x < 0) return float.NaN;
            if (float.IsPositiveInfinity(x)) return 0;
            return (float)TrigonometricIntegral((Complex)x, false).Real;
        }
        /// <summary>
        /// Returns the value of the integral cosine.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Ci(Complex32 z)
        {
            return (Complex32)TrigonometricIntegral((Complex)z, false);
        }
        /// <summary>
        /// Returns the value of the integral sine.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Si(float x)
        {
            if (float.IsInfinity(x)) return (float)(Math.Sign(x) * Math.PI / 2);
            return (float)TrigonometricIntegral((Complex)x, true).Real;
        }
        /// <summary>
        /// Returns the value of the integral sine.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Si(Complex32 z)
        {
            return (Complex32)TrigonometricIntegral((Complex)z, true);
        }
        /// <summary>
        /// Returns the value of an integral exponential function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Ei(float x)
        {
            if (float.IsNegativeInfinity(x)) return 0;
            if (float.IsPositiveInfinity(x)) return float.PositiveInfinity;
            return (float)ExponentialIntegral((Complex)x).Real;
        }
        /// <summary>
        /// Returns the value of an integral exponential function.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Ei(Complex32 z)
        {
            return (Complex32)ExponentialIntegral((Complex)z);
        }
        /// <summary>
        /// Returns the value of the integral logarithm.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Li(float x)
        {
            if (x < 0) return float.NaN;
            if (x == 0) return 0;
            return (float)ExponentialIntegral(Complex.Log((Complex)x)).Real;
        }
        /// <summary>
        /// Returns the value of the integral logarithm.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Li(Complex32 z)
        {
            if (z == 0) return Complex32.Zero;
            return (Complex32)ExponentialIntegral(Complex.Log((Complex)z));
        }
        #endregion

        #region Bessel functions

        /// <summary>
        /// Returns the value of a Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float J(float x, int a)
        {
            return (float)BesselJ((Complex)x, a).Real;
        }
        /// <summary>
        /// Returns the value of a Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 J(Complex32 x, int a)
        {
            return (Complex32)BesselJ((Complex)x, a);
        }

        /// <summary>
        /// Returns the value of a Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Y(float x, int a)
        {
            if (x < 0) return float.NaN;
            return (float)BesselY((Complex)x, a).Real;
        }
        /// <summary>
        /// Returns the value of a Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Y(Complex32 x, int a)
        {
            return (Complex32)BesselY((Complex)x, a);
        }

        /// <summary>
        /// Returns the value of the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float I(float x, int a)
        {
            return (float)BesselI((Complex)x, a).Real;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 I(Complex32 x, int a)
        {
            return (Complex32)BesselI((Complex)x, a);
        }

        /// <summary>
        /// Returns the value of the modified Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float K(float x, int a)
        {
            if (x < 0) return float.NaN;
            return (float)BesselK((Complex)x, a).Real;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 K(Complex32 x, int a)
        {
            return (Complex32)BesselK((Complex)x, a);
        }



        #endregion

        #region Struve functions
        /// <summary>
        /// Returns the value of the Struve function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float H(float x, int a)
        {
            return (float)StruveValue((Complex)x, a, false).Real;
        }
        /// <summary>
        /// Returns the value of the Struve function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 H(Complex32 x, int a)
        {
            return (Complex32)StruveValue((Complex)x, a, false);
        }
        /// <summary>
        /// Returns the value of the modified Struve function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="v">Value</param>
        /// <returns>Value</returns>
        public static float L(float x, int v)
        {
            return (float)StruveValue((Complex)x, v, true).Real;
        }
        /// <summary>
        /// Returns the value of the modified Struve function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="v">Value</param>
        /// <returns>Value</returns>
        public static Complex32 L(Complex32 x, int v)
        {
            return (Complex32)StruveValue((Complex)x, v, true);
        }
        #endregion

        #region Fibonacci & Lucas numbers
        /// <summary>
        /// Returns the value of the Fibonacci number.
        /// </summary>
        /// <param name="n">Integer number</param>
        /// <returns>Integer number</returns>
        public static int Fibonacci(int n)
        {
            return FibonacciValue(n, false);
        }
        /// <summary>
        /// Returns the value of the Luca number.
        /// </summary>
        /// <param name="n">Integer number</param>
        /// <returns>Integer number</returns>
        public static int Lucas(int n)
        {
            return FibonacciValue(n, true);
        }
        #endregion

        #region Harmonic number
        /// <summary>
        /// Returns the harmonic number.
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Value</returns>
        public static float Harm(int n)
        {
            if (n < 0) return float.NaN;
            if (n == 0) return 0;
            return (float)(Polygamma((Complex)((double)n + 1), false).Real + EulerGamma);
        }
        /// <summary>
        /// Returns the harmonic number.
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="m">Value</param>
        /// <returns>Value</returns>
        public static float Harm(int n, float m)
        {
            if (n < 0) return float.NaN;
            double sum = 0;
            for (int i = 1; i <= n; i++) sum += Math.Pow(i, -(double)m);
            return (float)sum;
        }
        #endregion

        #region Euler function
        /// <summary>
        /// Returns the Euler number.
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Value</returns>
        public static double Euler(int n)
        {
            // special cases:
            if (n < 0)
                return float.NaN;
            else if (n == 0)
                return 1;

            // for even number:
            else if (Maths.Mod(n, 2) == 0)
            {
                // get it from memory
                if (n <= 186)
                {
                    return Special.A122045[n / 2 - 1];
                }
                return float.NaN;
            }
            // for odd number:
            return 0;
        }
        /// <summary>
        /// Returns the value of the Euler polynomial.
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static double Euler(int n, float x)
        {
            return NumberPolynomial(n, x, true);
        }
        #endregion

        #region Bernoulli function
        /// <summary>
        /// Returns the Bernoulli number.
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Value</returns>
        public static double Bernoulli(int n)
        {
            // special cases:
            if (n < 0)
                return float.NaN;
            else if (n == 0)
                return 1;
            else if (n == 1)
                return -0.5f;

            // for even number:
            else if (Maths.Mod(n, 2) == 0)
            {
                // get it from memory
                if (n <= 258)
                {
                    return Special.A027641[n / 2 - 1];
                }
                return float.NaN;
            }
            // for odd number:
            return 0;
        }
        /// <summary>
        /// Returns the value of the Bernoulli polynomial.
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static double Bernoulli(int n, float x)
        {
            return NumberPolynomial(n, x, false);
        }
        #endregion

        #region Minkowski function
        /// <summary>
        /// Returns the value of the Minkowski function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Minkowski(long x)
        {
            // The question-mark function fixes every integer.
            return x;
        }
        #endregion

        #region Sequences of numbers
        /// <summary>
        /// Sequence A027641.
        /// </summary>
        private static readonly double[] A027641 = new double[]
        {
            0.166666666666666666666666666666666666666666666666666666666667,
           -0.0333333333333333333333333333333333333333333333333333333333333,
            0.0238095238095238095238095238095238095238095238095238095238095,
           -0.0333333333333333333333333333333333333333333333333333333333333,
            0.0757575757575757575757575757575757575757575757575757575757576,
           -0.253113553113553113553113553113553113553113553113553113553114,
            1.16666666666666666666666666666666666666666666666666666666667,
           -7.09215686274509803921568627450980392156862745098039215686275,
            54.9711779448621553884711779448621553884711779448621553884712,
           -529.1242424242424249314353801310062408447265625,
            6192.1231884057970091816969215869903564453125,
           -86580.253113553117145784199237823486328125,
            1425517.166666666977107524871826171875,
           -27298231.0678160898387432098388671875,
            601580873.90064239501953125,
           -15116315767.092159271240234375,
            429614643061.16668701171875,
           -13711655205088.330078125,
            488332318973593.1875,
           -19296579341940076.0,
            841693047573682560.0,
           -40338071854059462656.0,
            2115074863808199000064.0,
           -120866265222965295579136.0,
            7500866746076964166041600.0,
           -503877810148106884987486208.0,
            36528776484818117877519351808.0,
           -2849876930245088236122601947136.0,
            238654274996836310508426091298816.0,
           -21399949257225330247800505298321408.0,
            2050097572347810034157746982687342592.0,
           -209380059113463793012647415640108302336.0,
            22752696488463519698561496723250716082176.0,
           -2625771028623958030214791074126732688621568.0,
            321250821027180317428435351790225246396612608.0,
           -41598278166794711978200276325512922486305456128.0,
            5692069548203528317355790461742945664460956631040.0,
           -821836294197845776649876697189833540262035615383552.0,
            125029043271669897506885706274053657567199045890342912.0,
           -20015583233248370051792568366183437589101736793292668928.0,
            3367498291536436858725459209336816257030639921745753014272.0,
           -594709705031354502049750498053645919053879252111798906650624.0,
            1.1011910323627979045158206026592429963897952330719731013963e62,
           -2.13552595452534991558334720305307346520872020953392665238899e64,
            4.33288969866411863638232010122985676887721481700185983955699e66,
           -9.18855282416693318106533905537448294305383705378245339289825e68,
            2.0346896776329068803647719144955370450922522378940306576789e71,
           -4.7003833958035730157674974015948555777272519594093565148942e73,
            1.1318043445484249410725817313660083435449470330104860495433e76,
           -2.83822495706937114678248182940911866588389901514318497220466e78,
            7.40642489796788529354748598865042207929188648002583089441166e80,
           -2.0096454802756599870081060791981766846284745175520638572193e83,
            5.66571700508059420887233268977148771879706673802185161586432e85,
           -1.6584511154136219437234638454587009914837240155034651253857e88,
            5.03688599504923783896478298747871572448468662814872806462767e90,
           -1.58614682376581902614587000198040099996503141997053445150864e93,
            5.17567436175456311178752167302462403726113499716915409840838e95,
           -1.74889218402171188816021790530358849504622455237801259151532e98,
            6.11605199949521901296298005939495641146999201388699081095155e100,
           -2.21227769127078292149857803770811512559622750886224977635019e103,
            8.27227767987709687660613776431527490370915404345425054574528e105,
           -3.19589251114157084776409510523994112936024681625182396303699e108,
            1.27500822233877904696082368968356318782009773409266045580098e111,
           -5.25009230867741313467985686318449353023207998957580283537145e113,
            2.23018178942416301340661950936401895969403511914125731347591e116,
           -9.76845219309552074220408445655064473588790064695049483684759e118,
            4.40983619784529498205698125495662967577923441618079301678741e121,
           -2.05085708864640889574357632168670075648446415536760953855799e124,
            9.82144332797912765807607436077982175238879897214938800210545e126,
           -4.84126007982088805691451886595456745892784124244738914551734e129,
            2.45530888014809791682528532363557402092699052545481697315934e132,
           -1.28069268040847507813723639564508916744519481916574331612341e135,
            6.86761671046685794352718442702994340727704715936764246625598e137,
           -3.7846468581969101621095865768008278124904893666142088748761e140,
            2.14261012506652906510450120510522432260656756567484737427257e143,
           -1.24567271371836998349634545874581111163632446185107508042324e146,
            7.43457875510001443956855755955050576889972576151039708030967e148,
           -4.55357953046416979145968384968369983669462955623418794863421e151,
            2.86121128168588704603487881414828974143342817862557405123447e154,
           -1.84377235520338699519589434891927561212324166592681993388664e157,
            1.21811545362210496656724882070453175344651444513998349133864e160,
           -8.24821871853141216785484605290466498979217043473108315199036e162,
            5.72258779378329421678115626095616333281148998620922488296459e165,
           -4.06685305250590962105440857809913259792757881869118006357394e168,
            2.95960920646420479007743570819810233325436435218970900216399e171,
           -2.20495225651894615368307670360306991162801468142062178548309e174,
            1.68125970728895993749371058466736716170847214645904175801679e177,
           -1.31167362135569603944210519904485812542259573925565974960953e180,
            1.04678940094780397015084866777234879516265854029413770056151e183,
           -8.5432893578833710248124005359607211263683035276171302857864e185,
            7.1287821322486534488167835336906013750768485777448447528409e188,
           -6.08029314555359049665909285616042069926815131673966820380177e191,
            5.29967764248499211671611060091570164542757536816790348474136e194,
           -4.71942591687458604699269861534033695716707492189401616316378e197,
            4.29284137914029826070841811734553585864646091167302065201168e200,
           -3.98767449682322123296935998177114008533400284780121330523514e203,
            3.78197804193588797657318519554731981038162417744000041849395e206,
           -3.66142336836811918986123605424717958063923237299248880140341e209,
            3.61760902723728599167245528727089222102744767716779676157597e212,
           -3.64707726451913486561056842871572566300567152741905781213389e215,
            3.75087554364544063232692224124108557030929340399913341443501e218,
           -3.93458672964390300601992594022661357623121811698789216400172e221,
            4.20882111481900839569739650696605904148808088794013041370858e224,
           -4.59022962206179175425292693975421245971403882821905327048065e227,
            5.10317257726295805255711081303501141675201667482865380391774e230,
           -5.7822762303656963031228675308956212524549937115474674565998e233,
            6.67624821678358838017942792421516821882439233675200269178361e236,
           -7.85353076444504170220996386788735260526636756626106287201836e239,
            9.41068940670587242921381358271943326803708771306935461491816e242,
           -1.14849338734651807793649462205653591784510178822997246842936e246,
            1.42729587428487908796653191402977414634523039888257487145566e249,
           -1.80595595869093093539050272087441705873491055055840196644337e252,
            2.32615353076608109534439608243373828023729180025874119899477e255,
           -3.04957517154995898315383590877217052394096070579890754275451e258,
            4.06858060764339669246452724551233233088133189600849084296728e261,
           -5.52310313219743620065867776486536681636784118287796576464152e264,
            7.6277279396434385985638431216284788999441215837622777404542e267,
           -1.07155711196978895453618756016866974621961567085047422210229e271,
            1.53102008959691899726269085961622679829880291107894743575103e274,
           -2.22448916821798317978828962722646592764530144029166526492042e277,
            3.28626791906901399061342175480150778399861978247561459988121e280,
           -4.93559289559603397335492120298902934912343253048333017180432e283,
            7.53495712008325109008001025561176000306326320872542116983921e286,
           -1.16914851545841800139338970344528198850962357471005557092085e290,
            1.84352614678389384190513271871440592469742307916285568126251e293,
           -2.95368261729680817947377711654550539612299275825946981985138e296,
            4.80793212775015680793525997362458012705721533354406920469416e299,
           -7.95021250458852516728999753351542387618851025993657279302737e302,
            1.33527841873546293819631624429361412259570631983417101322206e306,
        };
        /// <summary>
        /// Sequence A000142.
        /// </summary>
        private static readonly double[] A000142 = new double[]
        {
            1.0,
            1.0,
            2.0,
            6.0,
            24.0,
            120.0,
            720.0,
            5040.0,
            40320.0,
            362880.0,
            3628800.0,
            39916800.0,
            479001600.0,
            6227020800.0,
            87178291200.0,
            1307674368000.0,
            20922789888000.0,
            355687428096000.0,
            6402373705728000.0,
            121645100408832000.0,
            2432902008176640000.0,
            51090942171709440000.0,
            1124000727777607680000.0,
            25852016738884978212864.0,
            620448401733239409999872.0,
            15511210043330986055303168.0,
            403291461126605650322784256.0,
            10888869450418351940239884288.0,
            304888344611713836734530715648.0,
            8841761993739700772720181510144.0,
            265252859812191032188804700045312.0,
            8222838654177922430198509928972288.0,
            263130836933693517766352317727113216.0,
            8683317618811885938715673895318323200.0,
            295232799039604119555149671006000381952.0,
            10333147966386144222209170348167175077888.0,
            371993326789901177492420297158468206329856.0,
            13763753091226343102992036262845720547033088.0,
            523022617466601037913697377988137380787257344.0,
            20397882081197441587828472941238084160318341120.0,
            815915283247897683795548521301193790359984930816.0,
            33452526613163802763987613764361857922667238129664.0,
            1.4050061177528797887796357975907848321789726105272e51,
            6.0415263063373834074440829285578945930237590418489e52,
            2.6582715747884485291342130280962418892431502625294e54,
            1.1962222086548018857499272315746937350318626585858e56,
            5.502622159812088456668950435842974564586819473163e57,
            2.5862324151116817767349100665299702655232519982624e59,
            1.2413915592536072528327568319343857274511609591659e61,
            6.0828186403426752248860160811673162316877754210242e62,
            3.0414093201713375576366966406747986832057064836515e64,
            1.551118753287382189470754582685817365323346291853e66,
            8.0658175170943876845634591553351679477960544579306e67,
            4.274883284060025484791254765342395718256495012315e69,
            2.3084369733924137924371883906026708550254478496563e71,
            1.2696403353658276446882823840816011312245221598828e73,
            7.1099858780486348102543813508569663348573240953439e74,
            4.0526919504877220527556156789809444757511993541236e76,
            2.3505613312828789062977962804562476349569662739554e78,
            1.3868311854568986493322118514385335285353380986813e80,
            8.3209871127413915800563961029596410774579455410767e81,
            5.0758021387722483583354016137308849072428138984387e83,
            3.1469973260387939390320343330721249710233204778006e85,
            1.982608315404440084965732774767545707658109829136e87,
            1.2688693218588416543780689758512292529011902906471e89,
            8.2476505920824715167353803272950205238422572101466e90,
            5.4434493907744306944549606027563585676128303456872e92,
            3.6471110918188683221214362054827498508015278133658e94,
            2.4800355424368305479709011539871079838475553997611e96,
            1.7112245242814129737573543427207344887665272148063e98,
            1.1978571669969890269925854460558840225267029209529e100,
            8.5047858856786217613936449886228345036310665787676e101,
            6.1234458376886076682034243918084408426143679367127e103,
            4.4701154615126833670305181118791598550111254536754e105,
            3.3078854415193855897507845860662739792859484152509e107,
            2.4809140811395391401649674453868616759922516881581e109,
            1.8854947016660498466497675672866749860207537596979e111,
            1.4518309202828583792503372319096021362422032622554e113,
            1.132428117820629460628519376473454765964154487391e115,
            8.9461821307829729139453610567812466009169986197986e116,
            7.1569457046263778832073404098641551692451427821501e118,
            5.7971260207473655478592076093169551533024183171149e120,
            4.7536433370128398180495087193420485740326098790968e122,
            3.9455239697206569095363763848575524105091557652835e124,
            3.3142401345653519918939627851870022559861385859851e126,
            2.817104114380549361453400700637317692706186975947e128,
            2.4227095383672724277628115968482030522825707406365e130,
            2.1077572983795269087233798237224287232533562814768e132,
            1.85482642257398355359015441641340379717002520724e134,
            1.6507955160908452497218052643056785820348586593118e136,
            1.4857159644817606885981264446583904686485043853385e138,
            1.3520015276784022811614124898346447400031519033616e140,
            1.2438414054641300055918190849808704283732243800785e142,
            1.1567725070816408727081687710539103627025574530463e144,
            1.0873661566567424099481672918376200175213508153269e146,
            1.0329978488239052206885505130495304991006115078121e148,
            9.916779348709491027158490294784105978020630335002e149,
            9.6192759682482062236598631563798937437476306361515e151,
            9.4268904488832420294101483608740343761375924661801e153,
            9.3326215443944096091160468772652940323762165415183e155,
            9.3326215443944102188325606108575267240944254854961e157,
            9.4259477598383536381383908353428013766109757230959e159,
            9.6144667150351210855109846899687251699348628127378e161,
            9.9029007164861753574104173353829959119840213587604e163,
            1.0299016745145621553359182054762848245165958006211e166,
            1.0813967582402903482108699210497876860853577081697e168,
            1.1462805637347078319526217879186988515037213497498e170,
            1.2265202031961373185133888353370611130684668524854e172,
            1.3246418194518283589128412232089037219691505933354e174,
            1.4438595832024928189521163811423104875896274070817e176,
            1.5882455415227421289655392351515189289144558208272e178,
            1.762952551090243665975210588885144142387414195101e180,
            1.974506857221072832182032249755631903506040988586e182,
            2.2311927486598122561395742763464263293811085711459e184,
            2.543559733472186205984785013970489854470847792379e186,
            2.9250936934930141417131746698376362635106648249008e188,
            3.3931086844518965033194432062534716898733796047438e190,
            3.9699371608087190355169141055460833161291448996096e192,
            4.6845258497542883302114664681456760164960085889976e194,
            5.5745857612076033333946596938961053806019179795001e196,
            6.6895029134491239336813425305794391205237757563415e198,
            8.0942985252734400128682237436778312047925893204201e200,
            9.8750442008335975771386594694042667512406118845316e202,
            1.2146304367025324845837253661169005205421689094861e205,
            1.5061417415111403610829709356251069731491340792132e207,
            1.8826771768889253829171044051984554368188435135121e209,
            2.3721732428800459167764066567904427019588717049411e211,
            3.0126600184576582348830701825566878653822517532445e213,
            3.8562048236258025406503298336725604676892822441529e215,
            4.9745042224772854994195754995984321775186454084691e217,
            6.4668554892204716391337792152122055255592791711738e219,
            8.4715806908788174208664474123129035265417567762389e221,
            1.118248651196003918817564057259678902372371328839e224,
            1.4872707060906851873704731576473721249631137913553e226,
            1.9929427461615181195156186219572376040046331831353e228,
            2.6904727073180495455082595644162468812620929233147e230,
            3.659042881952547209527099785668992673243689893492e232,
            5.0128887482749898425216945994261389242057013770114e234,
            6.9177864726194858697193801988793718374395142120906e236,
            9.6157231969410858830469292126874942877955258679849e238,
            1.3462012475717519819847898602483572772785942779053e241,
            1.8981437590761701317329243607966509065128364990507e243,
            2.6953641378881613975443839032798606020700676159442e245,
            3.8543707171800705947204693672221529226566388300967e247,
            5.5502938327393012589528442598184128104868108689386e249,
            8.0479260574719866188104157296663251281737262558327e251,
            1.1749972043909099283211628879889409709620811165324e254,
            1.7272458904546375894227231696563570756893750298017e256,
            2.5563239178728636856897457541302426646590893698766e258,
            3.8089226376305670581113614183352242513751438575187e260,
            5.7133839564458504888431685060296571931908532210186e262,
            8.6272097742332346157168591505617904277861406143342e264,
            1.3113358856834517618251263055477563604892371412407e267,
            2.0063439050956811222768492161921665482364459558907e269,
            3.0897696138473488989801101804175762109592920239887e271,
            4.7891429014633911950713211206987054304438393469277e273,
            7.471062926282890533380664092650436713658681188259e275,
            1.1729568794264138444915531933301707060062748673626e278,
            1.8532718694937337798302304500930615353830744697171e280,
            2.9467022724950368550304930581341500059031943714696e282,
            4.7147236359920589680487888930146400094451109943514e284,
            7.5907050539472147668628274986623653506961986030378e286,
            1.2296942187394487685641830599055001502279618094322e289,
            2.0044015765453015187263511820261834335817828575911e291,
            3.2872185855342944907112159385229408310741238864494e293,
            5.4239106661315859522174013088754018516551044619436e295,
            9.0036917057784328985056286255334204133074096592532e297,
            1.5036165148649982533898213892747005856378159793089e300,
            2.5260757449731968991490061844697939504885310524502e302,
            4.2690680090047026720062648805820850387756975875746e304,
            7.2574156153079940453996357155895914678961841172423e306,
        };
        /// <summary>
        /// Sequence A122045.
        /// </summary>
        private static readonly double[] A122045 = new double[]
        {
            -1.0,
             5.0,
            -61.0,
             1385.0,
            -50521.0,
             2702765.0,
            -199360981.0,
             19391512145.0,
            -2404879675441.0,
             370371188237525.0,
            -69348874393137904.0,
             15514534163557089280.0,
            -4087072509293124124672.0,
             1252259641403629925040128.0,
            -441543893249023112372027392.0,
             177519391579539304507368275968.0,
            -80723299235887898068046850293760.0,
             41222060339517699219515483317338112.0,
            -23489580527043111237325070972959588352.0,
             14851150718114980007771290846864380788736.0,
            -10364622733519610322736028694318548613332992.0,
             7947579422597593581036447205938820399155380224.0,
            -6667537516685544830982458769181946136668143616000.0,
             6096278645568541790228442449457750220561157371985920.0,
            -6053285248188621883417017917245308857755712373867413504.0,
             6506162486684608510440254565853912312843020577732143611904.0,
            -7546659939008739271231891941068505271100907148104580850515977.0,
             9.42032189642024117029653987727040357142967393494904530837648e63,
            -1.26220192518062200766940816422423442090824669296392239118522e67,
             1.81089114965792304567172651019370789035178011401985350135716e70,
            -2.77571017020715785153706725905409592784929360042138526804035e73,
             4.53581033300178873183845375182440770271317432262387897678718e76,
            -7.88628420666178909321870785703324111843713756139283816015288e79,
             1.45618443801396308903980797707391570583250082893301288469e83,
            -2.85051783223697701075758089194738584661800764710150138510653e86,
             5.90574720777544426844910678644178262854134720057761800475874e89,
            -1.29297366418786410389269648309290712538591911532289247513896e93,
             2.98692818328457691000610068544311841203205288807302728536229e96,
            -7.27060171401686442513144180750678121046414109797425110784264e99,
             1.8622915758412701384349147366598671113782545828754552066129e103,
            -5.01310494081097996484468299887663804251446455532450853954811e106,
             1.41652557597856293895025676203809322008072905892279836863327e110,
            -4.1966431640402452779591075681258768351693864961916057933188e113,
             1.302159590524046126369484063896183128865886094972353427233e117,
            -4.22724068613990882320942470716279737275467600092365310804626e120,
             1.43432127919765801478114221331805669282627309913289094288186e124,
            -5.08179907245804290327013129412064438468021872302714025860577e127,
             1.87833293645293009429932609249457210582511142027270809253815e131,
            -7.23653438103385801091977060254655023485376375699570171559289e134,
             2.90352834666109691193648582021635560320923721009836362684885e138,
            -1.21229373789292196985790359148899045066488864191107837624119e142,
             5.26306424961699132922867650270045256917210200404591849896002e145,
            -2.37407307193676613444397971185801838736518087157756429507812e149,
             1.11189009424828205369640983214229021354876102213780649821977e153,
            -5.40307865979529323090987512189247587858209495258119340185774e156,
             2.72234108557222712607694470521316337763813999968686359468173e160,
            -1.42130105480096708921724268915539592541828996834320545228379e164,
             7.68426182064690282266660093248723270977854299159375532739554e167,
            -4.29962192543975040028581899870282367007480315711551073982038e171,
             2.48839157478298713619604697758333796024171557967378371972775e175,
            -1.48875820890620386196896135481729129349058104526270299753085e179,
             9.20261411885209372269704777166644446959768272803345768890992e182,
            -5.87424445729243577521888833616725549721731744854197071892569e186,
             3.87013355417592696733856268773112115587056868765773113470594e190,
            -2.63038464627282205142925192000387504737523085666729635567083e194,
             1.84342186190681611237463552747621795310801031741771707310489e198,
            -1.33150076083199794467135216187391299040926674495401811176873e202,
             9.90773407946409988676033835119472855742395813525809239244951e205,
            -7.59161615376086661069934157407387899531523434692866129479005e209,
             5.98738690421595488843273439931836679779298506469340158923842e213,
            -4.85853153680526962878834959532910209742077959025966808866591e217,
             4.05474737750791489372344949537425177494580382647116914394292e221,
            -3.47892371339090604367719331387826727192214732827580244097765e225,
             3.06749738825108493077753427507150071973348947259441411043398e229,
            -2.77857404780457398412720199541452677766447630875005472975504e233,
             2.58465603902711788716741587393145364906078903144675176202619e237,
            -2.46817048046364094667210893288509701822698969438155373850515e241,
             2.4187539760367128249917153087475578880169651392238054590364e245,
            -2.43169264709107282615088875463734917330677608026924951798079e249,
             2.50718300057371384087291375129517485752197239767889561857948e253,
            -2.6502520005258141935263902742471414084327859710059871475327e257,
             2.87130197316668000171507992791938863839734029604874444880572e261,
            -3.18736021623541101072259543193584569386177442843757875358645e265,
             3.62424164505845579420550817104847545918510713215729565950763e269,
            -4.22000551313026087233885754179495565508060473653385418256369e273,
             5.03034557853149999562618903891101674615476791956764118621342e277,
            -6.13696178494213373370421222447301459117641511081674038800514e281,
             7.66062813846337350526483265240843883880613843320540624512461e285,
            -9.78178011283967467592520364785782561305562032336418074186079e289,
             1.27733166367198099806106221718741466821274419712357185992648e294,
            -1.70535141854472089452708066596252042944127306749367113463532e298,
             2.32725003482002998572779903012518091723702509768829825693762e302,
            -3.245547458389247260235183311044519483521857479024361599136e306
        };
        #endregion
    }
}
