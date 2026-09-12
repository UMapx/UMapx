using System;
using System.Numerics;

namespace UMapx.Core
{
    /// <summary>
    /// Used to implement special mathematical functions.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// <see href="https://en.wikipedia.org/wiki/Special_functions"/>.
    /// </remarks>
    public static partial class Special
    {
        #region Chebyshev polynomial
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
        public static float ChebyshevT(float x, int n)
        {
            return (float)ChebyshevValue((Complex)x, n, false).Real;
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
        public static Complex32 ChebyshevT(Complex32 x, int n)
        {
            return (Complex32)ChebyshevValue((Complex)x, n, false);
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the second kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
        public static float ChebyshevU(float x, int n)
        {
            return (float)ChebyshevValue((Complex)x, n, true).Real;
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the second kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
        public static Complex32 ChebyshevU(Complex32 x, int n)
        {
            return (Complex32)ChebyshevValue((Complex)x, n, true);
        }
        #endregion

        #region Abel polynomial
        /// <summary>
        /// Returns the value of the Abel polynomial.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Power.</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
        public static float Abel(float x, float a, int n)
        {
            if (n < 0) return float.NaN;
            if (n == 0) return 1;
            return (float)((Complex)x * Complex.Pow((Complex)x - n * (Complex)a, n - 1)).Real;
        }
        /// <summary>
        /// Returns the value of the Abel polynomial.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Complex power.</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
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
        /// <param name="x">Value.</param>
        /// <param name="a">Power.</param>
        /// <param name="k">Order.</param>
        /// <returns>Value.</returns>
        public static float Laguerre(float x, float a, int k)
        {
            return (float)OrthogonalPolynomial((Complex)x, (Complex)a, k, 0).Real;
        }
        /// <summary>
        /// Returns the value of the Laguerre polynomial.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Power.</param>
        /// <param name="k">Order.</param>
        /// <returns>Value.</returns>
        public static Complex32 Laguerre(Complex32 x, Complex32 a, int k)
        {
            return (Complex32)OrthogonalPolynomial((Complex)x, (Complex)a, k, 0);
        }
        #endregion

        #region Legendre polynomial
        /// <summary>
        /// Returns the value of the Legendre polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="m">Order.</param>
        /// <returns>Value.</returns>
        public static float Legendre(float x, int m)
        {
            return (float)OrthogonalPolynomial((Complex)x, Complex.Zero, m, 1).Real;
        }
        /// <summary>
        /// Returns the value of the Legendre polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="m">Order.</param>
        /// <returns>Value.</returns>
        public static Complex32 Legendre(Complex32 x, int m)
        {
            return (Complex32)OrthogonalPolynomial((Complex)x, Complex.Zero, m, 1);
        }
        #endregion

        #region Hermite polynomial
        /// <summary>
        /// Returns the value of the Hermite polynomial.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="m">Order.</param>
        /// <returns>Value.</returns>
        public static float Hermite(float x, int m)
        {
            return (float)OrthogonalPolynomial((Complex)x, Complex.Zero, m, 2).Real;
        }
        /// <summary>
        /// Returns the value of the Hermite polynomial.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="m">Order.</param>
        /// <returns>Value.</returns>
        public static Complex32 Hermite(Complex32 x, int m)
        {
            return (Complex32)OrthogonalPolynomial((Complex)x, Complex.Zero, m, 2);
        }
        #endregion

        #region Gegenbauer polynomial
        /// <summary>
        /// Returns the value of the Gegenbauer polynomial.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Power.</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
        public static float Gegenbauer(float x, float a, int n)
        {
            return (float)OrthogonalPolynomial((Complex)x, (Complex)a, n, 3).Real;
        }
        /// <summary>
        /// Returns the value of the Gegenbauer polynomial.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Power.</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
        public static Complex32 Gegenbauer(Complex32 x, Complex32 a, int n)
        {
            return (Complex32)OrthogonalPolynomial((Complex)x, (Complex)a, n, 3);
        }
        #endregion

        #region Sinc function
        /// <summary>
        /// Returns the value of the normalized cardinal sine function: f(x) = sin(πx) / (πx).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Sinc(float x)
        {
            return Special.Sinc(x, Maths.Pi);
        }
        /// <summary>
        /// Returns the value of the normalized cardinal sine function: f(x) = sin(πx) / (πx).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Sinc(Complex32 x)
        {
            return Special.Sinc(x, Maths.Pi);
        }
        /// <summary>
        /// Returns the value of the cardinal sine function with the parameter: f(x, a) = sin(ax) / (ax).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
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
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
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
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Agd(float x)
        {
            // gd^{-1}(x) = artanh(sin(x))
            return Maths.Atanh(Maths.Sin(x));
        }
        /// <summary>
        /// Returns the value of the inverse Guderman function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Agd(Complex32 x)
        {
            // gd^{-1}(x) = artanh(sin(x))
            return Maths.Atanh(Maths.Sin(x));
        }
        /// <summary>
        /// Returns the value of the Guderman function.
        /// </summary>
        /// <param name="x">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Gd(float x)
        {
            return (float)(2 * Math.Atan(Math.Tanh((double)x / 2)));
        }
        /// <summary>
        /// Returns the value of the Guderman function.
        /// </summary>
        /// <param name="x">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static Complex32 Gd(Complex32 x)
        {
            return (Complex32)(2 * Complex.Atan(2 * LogisticValue((Complex)x) - 1));
        }
        /// <summary>
        /// Returns the value of the function Cas(x).
        /// </summary>
        /// <param name="theta">Theta.</param>
        /// <returns>Value.</returns>
        public static float Cas(float theta)
        {
            return Maths.Cos(theta) + Maths.Sin(theta);
        }
        /// <summary>
        /// Returns the value of the function Cas(x).
        /// </summary>
        /// <param name="theta">Theta.</param>
        /// <returns>Value.</returns>
        public static Complex32 Cas(Complex32 theta)
        {
            return Maths.Cos(theta) + Maths.Sin(theta);
        }
        #endregion

        #region Rademacher function
        /// <summary>
        /// Returns the value of the Radamecher function.
        /// </summary>
        /// <param name="t">Value [0, 1].</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
        public static float Rademacher(float t, int n)
        {
            return (float)RademacherValue((Complex)t, n).Real;
        }
        /// <summary>
        /// Returns the value of the Radamecher function.
        /// </summary>
        /// <param name="z">Value.</param>
        /// <param name="n">Order.</param>
        /// <returns>Value.</returns>
        public static Complex32 Rademacher(Complex32 z, int n)
        {
            return (Complex32)RademacherValue((Complex)z, n);
        }
        #endregion

        #region Heavyside delta-function
        /// <summary>
        /// Returns the value of the Heaviside delta function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="k">Smoothing factor.</param>
        /// <returns>Value.</returns>
        public static float Heaviside(float x, float k)
        {
            return (float)LogisticValue(2 * (Complex)k * (Complex)x).Real;
        }
        /// <summary>
        /// Returns the value of the Heaviside delta function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="k">Smoothing factor.</param>
        /// <returns>Value.</returns>
        public static Complex32 Heaviside(Complex32 x, Complex32 k)
        {
            return (Complex32)LogisticValue(2 * (Complex)k * (Complex)x);
        }
        #endregion

        #region Mahler function
        /// <summary>
        /// Returns the value of the Mahler function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="t">Value.</param>
        /// <returns>Value.</returns>
        public static float Mahler(float x, float t)
        {
            return Maths.Exp(x * (1.0f + t - Maths.Pow(Maths.E, t)));
        }
        /// <summary>
        /// Returns the value of the Mahler function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="t">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Mahler(Complex32 x, Complex32 t)
        {
            return Maths.Exp(x * (1.0f + t - Maths.Pow(Maths.E, t)));
        }
        #endregion

        #region Gompertz function
        /// <summary>
        /// Gets the value of the Gompertz function.
        /// </summary>
        /// <param name="t">Value.</param>
        /// <param name="a">Upper asymptote.</param>
        /// <param name="b">Growth parameter.</param>
        /// <param name="c">Growth rate.</param>
        /// <returns>Value.</returns>
        public static float Gompertz(float t, float a, float b, float c)
        {
            return a * Maths.Exp(-b * Maths.Exp(-c * t));
        }
        /// <summary>
        /// Gets the value of the Gompertz function.
        /// </summary>
        /// <param name="t">Value.</param>
        /// <param name="a">Upper asymptote.</param>
        /// <param name="b">Growth parameter.</param>
        /// <param name="c">Growth rate.</param>
        /// <returns>Value.</returns>
        public static Complex32 Gompertz(Complex32 t, Complex32 a, Complex32 b, Complex32 c)
        {
            return a * Maths.Exp(-b * Maths.Exp(-c * t));
        }
        #endregion

        #region Dirac delta-function
        /// <summary>
        /// Returns the value of the Dirac delta function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Coefficient.</param>
        /// <returns>Value.</returns>
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
        /// <param name="x">Value.</param>
        /// <param name="a">Coefficient.</param>
        /// <returns>Value.</returns>
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
        /// <param name="x">Value.</param>
        /// <param name="a">Lower asymptote.</param>
        /// <param name="k">Upper asymptote.</param>
        /// <param name="b">Growth rate.</param>
        /// <param name="v">Affect.</param>
        /// <param name="q">Central moment.</param>
        /// <param name="c">Offset.</param>
        /// <returns>Value.</returns>
        public static float Logistic(float x, float a, float k, float b, float v, float q, float c)
        {
            return a + (k - a) / Maths.Pow(c + q * Maths.Exp(-b * x), 1.0f / v);
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Lower asymptote.</param>
        /// <param name="k">Upper asymptote.</param>
        /// <param name="b">Growth rate.</param>
        /// <param name="v">Affect.</param>
        /// <param name="q">Central moment.</param>
        /// <param name="c">Offset.</param>
        /// <returns>Value.</returns>
        public static Complex32 Logistic(Complex32 x, Complex32 a, Complex32 k, Complex32 b, Complex32 v, Complex32 q, Complex32 c)
        {
            return a + (k - a) / Maths.Pow(c + q * Maths.Exp(-b * x), 1.0f / v);
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Lower asymptote.</param>
        /// <param name="k">Upper asymptote.</param>
        /// <param name="b">Growth rate.</param>
        /// <returns>Value.</returns>
        public static float Logistic(float x, float a, float k, float b)
        {
            return (float)((Complex)a + ((Complex)k - (Complex)a) * LogisticValue((Complex)b * (Complex)x)).Real;
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Lower asymptote.</param>
        /// <param name="k">Upper asymptote.</param>
        /// <param name="b">Growth rate.</param>
        /// <returns>Value.</returns>
        public static Complex32 Logistic(Complex32 x, Complex32 a, Complex32 k, Complex32 b)
        {
            return (Complex32)((Complex)a + ((Complex)k - (Complex)a) * LogisticValue((Complex)b * (Complex)x));
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Logistic(float x)
        {
            return (float)LogisticValue((Complex)x).Real;
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Logistic(Complex32 x)
        {
            return (Complex32)LogisticValue((Complex)x);
        }
        #endregion

        #region Elrang B and C functions
        /// <summary>
        /// Returns the value of the Erlang C-function.
        /// </summary>
        /// <param name="y">First parameter.</param>
        /// <param name="v">Second parameter.</param>
        /// <param name="t">Time parameter.</param>
        /// <returns>Value.</returns>
        public static float Erlang(float y, int v, float t)
        {
            Complex traffic = (Complex)y, blocking = ErlangBlocking(traffic, v);
            return (float)(v * blocking / (v - traffic + traffic * blocking) * Complex.Exp(-(v - traffic) * (Complex)t)).Real;
        }
        /// <summary>
        /// Returns the value of the Erlang C-function.
        /// </summary>
        /// <param name="y">First parameter.</param>
        /// <param name="v">Second parameter.</param>
        /// <param name="t">Time parameter.</param>
        /// <returns>Value.</returns>
        public static Complex32 Erlang(Complex32 y, int v, Complex32 t)
        {
            Complex traffic = (Complex)y, blocking = ErlangBlocking(traffic, v);
            return (Complex32)(v * blocking / (v - traffic + traffic * blocking) * Complex.Exp(-(v - traffic) * (Complex)t));
        }
        /// <summary>
        /// Returns the value of the Erlang B-function.
        /// </summary>
        /// <param name="y">First parameter.</param>
        /// <param name="v">Second parameter.</param>
        /// <returns>Value.</returns>
        public static float Erlang(float y, int v)
        {
            return (float)ErlangBlocking((Complex)y, v).Real;
        }
        /// <summary>
        /// Returns the value of the Erlang B-function.
        /// </summary>
        /// <param name="y">First parameter.</param>
        /// <param name="v">Second parameter.</param>
        /// <returns>Value.</returns>
        public static Complex32 Erlang(Complex32 y, int v)
        {
            return (Complex32)ErlangBlocking((Complex)y, v);
        }
        #endregion

        #region Lambert W-function
        /// <summary>
        /// Returns the value of the Lambert W-function.
        /// </summary>
        /// <param name="x">Value [-1/e,+inf).</param>
        /// <param name="k">Branch.</param>
        /// <returns>Value.</returns>
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
        /// <param name="z">Value.</param>
        /// <param name="k">Branch.</param>
        /// <returns>Value.</returns>
        public static Complex32 LambertW(Complex32 z, int k = 0)
        {
            return (Complex32)LambertValue((Complex)z, k);
        }
        /// <summary>
        /// Returns the value of the square super-root.
        /// </summary>
        /// <param name="x">Value [1,+inf).</param>
        /// <param name="k">Branch.</param>
        /// <returns>Value.</returns>
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
        /// <param name="z">Value.</param>
        /// <param name="k">Branch.</param>
        /// <returns>Value.</returns>
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
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Fresnelc(float x)
        {
            return (float)FresnelValue((Complex)x, false).Real;
        }
        /// <summary>
        /// Returns the value of the Fresnel integral C(x).
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Fresnelc(Complex32 z)
        {
            return (Complex32)FresnelValue((Complex)z, false);
        }
        /// <summary>
        /// Returns the value of the Fresnel integral S(x).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Fresnels(float x)
        {
            return (float)FresnelValue((Complex)x, true).Real;
        }
        /// <summary>
        /// Returns the value of the Fresnel integral S(x).
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Fresnels(Complex32 z)
        {
            return (Complex32)FresnelValue((Complex)z, true);
        }
        #endregion

        #region Owen's T-function
        /// <summary>
        /// Returns the value of the Owen T function.
        /// </summary>
        /// <param name="h">First value.</param>
        /// <param name="a">Second value.</param>
        /// <returns>Value.</returns>
        public static float Owen(float h, float a)
        {
            return (float)OwenValue((Complex)h, (Complex)a).Real;
        }
        /// <summary>
        /// Returns the value of the Owen T function.
        /// </summary>
        /// <param name="h">First value.</param>
        /// <param name="a">Second value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Owen(Complex32 h, Complex32 a)
        {
            return (Complex32)OwenValue((Complex)h, (Complex)a);
        }

        #endregion

        #region Riemann's Zeta function
        /// <summary>
        /// Returns the value of the Riemann zeta ζ(s) on the principal branch (real s).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <returns>Value.</returns>
        public static float Zeta(float s)
        {
            if (s == 1) return float.PositiveInfinity;
            if (float.IsPositiveInfinity(s)) return 1;
            return (float)ZetaValue((Complex)s).Real;
        }
        /// <summary>
        /// Returns the value of the Riemann zeta ζ(s) on the principal branch (complex s).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Zeta(Complex32 s)
        {
            return (Complex32)ZetaValue((Complex)s);
        }



        #endregion



        #region Gamma functions

        /// <summary>
        /// Returns the value of the Euler Gamma function: Gamma(z).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Gamma(float x)
        {
            return (float)GammaValue((double)x);
        }
        /// <summary>
        /// Returns the value of the Euler Gamma function: Gamma(z).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Gamma(Complex32 x)
        {
            return (Complex32)GammaValue((Complex)x);
        }

        /// <summary>
        /// Returns log(abs(Gamma(x))) for real arguments. Gamma poles return NaN.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float LogGamma(float x)
        {
            return (float)GammaLog((double)x);
        }
        /// <summary>
        /// Returns analytic log-gamma with its cut on the negative real axis.
        /// Its imaginary part is not reduced modulo 2*pi; the upper side is used on the cut.
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 LogGamma(Complex32 z)
        {
            return (Complex32)GammaLog((Complex)z);
        }

        /// <summary>
        /// Returns the value of the Digamma function: ψ(z).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float DiGamma(float x)
        {
            if (float.IsPositiveInfinity(x)) return float.PositiveInfinity;
            return (float)Polygamma((Complex)x, false).Real;
        }
        /// <summary>
        /// Returns the value of the Digamma function: ψ(z).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 DiGamma(Complex32 x)
        {
            return (Complex32)Polygamma((Complex)x, false);
        }

        /// <summary>
        /// Returns the value of the Trigamma function: ψ1(z).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float TriGamma(float x)
        {
            if (float.IsPositiveInfinity(x)) return 0;
            return (float)Polygamma((Complex)x, true).Real;
        }
        /// <summary>
        /// Returns the value of the Trigamma function: ψ1(z).
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 TriGamma(Complex32 x)
        {
            return (Complex32)Polygamma((Complex)x, true);
        }

        /// <summary>
        /// Returns the value of the incomplete upper Gamma function: Q(s, x) = Γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float GammaQ(float s, float x)
        {
            return (float)IncompleteGamma((double)s, (double)x, true, true);
        }
        /// <summary>
        /// Returns the value of the incomplete upper Gamma function: Q(s, x) = Γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 GammaQ(Complex32 s, Complex32 x)
        {
            return (Complex32)IncompleteGamma((Complex)s, (Complex)x, true, true);
        }

        /// <summary>
        /// Returns the value of an incomplete lower Gamma function: P(s, x) = γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float GammaP(float s, float x)
        {
            return (float)IncompleteGamma((double)s, (double)x, false, true);
        }
        /// <summary>
        /// Returns the value of an incomplete lower Gamma function: P(s, x) = γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 GammaP(Complex32 s, Complex32 x)
        {
            return (Complex32)IncompleteGamma((Complex)s, (Complex)x, false, true);
        }

        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float GammaIncomplete(float s, float x)
        {
            return (float)IncompleteGamma((double)s, (double)x, false, false);
        }
        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 GammaIncomplete(Complex32 s, Complex32 x)
        {
            return (Complex32)IncompleteGamma((Complex)s, (Complex)x, false, false);
        }

        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x) (complemented).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float GammaIncompleteComplemented(float s, float x)
        {
            return (float)IncompleteGamma((double)s, (double)x, true, false);
        }
        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x) (complemented).
        /// </summary>
        /// <param name="s">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 GammaIncompleteComplemented(Complex32 s, Complex32 x)
        {
            return (Complex32)IncompleteGamma((Complex)s, (Complex)x, true, false);
        }



        #endregion

        #region Generalized error function
        /// <summary>
        /// Returns n!/sqrt(pi) times the integral of exp(-t^n) from zero to x.
        /// </summary>
        /// <param name="x">Real integration endpoint.</param>
        /// <param name="n">Order [0, +inf).</param>
        /// <returns>Value.</returns>
        public static float Gerf(float x, int n)
        {
            return (float)GeneralizedErf((Complex)x, n).Real;
        }
        /// <summary>
        /// Returns the entire continuation of n!/sqrt(pi) times the integral of exp(-t^n).
        /// </summary>
        /// <param name="x">Complex integration endpoint.</param>
        /// <param name="n">Order [0, +inf).</param>
        /// <returns>Value.</returns>
        public static Complex32 Gerf(Complex32 x, int n)
        {
            return (Complex32)GeneralizedErf((Complex)x, n);
        }
        /// <summary>
        /// Returns the value of the generalized error function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Gerf(float x)
        {
            return Gerf(x, 2);
        }
        /// <summary>
        /// Returns the value of the generalized error function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Gerf(Complex32 x) => Gerf(x, 2);
        #endregion

        #region Factorial function
        /// <summary>
        /// Returns the natural logarithm of the factorial of a number log(n!).
        /// </summary>
        /// <param name="n">Value.</param>
        /// <returns>Value.</returns>
        public static float LogFactorial(float n)
        {
            return (float)GammaLog((double)n + 1);
        }
        /// <summary>
        /// Returns the natural logarithm of the factorial of a number log(n!).
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 LogFactorial(Complex32 z)
        {
            return (Complex32)GammaLog((Complex)z + 1);
        }
        /// <summary>
        /// Returns the factorial of a number.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <returns>Value.</returns>
        public static double Factorial(float n)
        {
            if (n >= 0 && n <= 170 && n == Math.Floor(n)) return A000142[(int)n];
            return GammaValue((double)n + 1);
        }
        /// <summary>
        /// Returns the factorial of a number.
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Factorial(Complex32 z)
        {
            return (Complex32)GammaValue((Complex)z + 1);
        }
        /// <summary>
        /// Returns the decreasing factorial of a number.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <param name="k">Value.</param>
        /// <returns>Value.</returns>
        public static double FactorialDown(float n, float k)
        {
            return FactorialProduct((Complex)n, (Complex)k, false).Real;
        }
        /// <summary>
        /// Returns the decreasing factorial of a number.
        /// </summary>
        /// <param name="z">Value.</param>
        /// <param name="k">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 FactorialDown(Complex32 z, Complex32 k)
        {
            return (Complex32)FactorialProduct((Complex)z, (Complex)k, false);
        }
        /// <summary>
        /// Returns the increasing factorial of a number (Pochhammer symbol).
        /// </summary>
        /// <param name="n">Value.</param>
        /// <param name="k">Value.</param>
        /// <returns>Value.</returns>
        public static float FactorialUp(float n, float k)
        {
            return (float)FactorialProduct((Complex)n, (Complex)k, true).Real;
        }
        /// <summary>
        /// Returns the increasing factorial of a number (Pochhammer symbol).
        /// </summary>
        /// <param name="z">Value.</param>
        /// <param name="k">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 FactorialUp(Complex32 z, Complex32 k)
        {
            return (Complex32)FactorialProduct((Complex)z, (Complex)k, true);
        }
        #endregion

        #region Binomial function
        /// <summary>
        /// Returns the value of binomial coefficients: C(n, k) = n! / k! / (n-k)! for k > 0.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <param name="k">Value.</param>
        /// <returns>Value.</returns>
        public static double Binomial(float n, float k)
        {
            return BinomialValue((Complex)n, (Complex)k).Real;
        }
        /// <summary>
        /// Returns the value of binomial coefficients: C(n, k) = n! / k! / (n-k)! for k > 0.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <param name="k">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Binomial(Complex32 n, Complex32 k)
        {
            return (Complex32)BinomialValue((Complex)n, (Complex)k);
        }
        /// <summary>
        /// Returns the natural logarithm of binomial coefficients: log(C(n, k)) = log(n!) - log(k!) - log(n-k!).
        /// </summary>
        /// <param name="n">Value.</param>
        /// <param name="k">Value.</param>
        /// <returns>Value.</returns>
        public static float LogBinomial(float n, float k)
        {
            if (k < 0 || (n >= 0 && n == Math.Floor(n) && k > n)) return float.NegativeInfinity;
            return (float)(GammaLog((double)n + 1) - GammaLog((double)k + 1) - GammaLog((double)n - k + 1));
        }
        /// <summary>
        /// Returns the natural logarithm of binomial coefficients: log(C(n, k)) = log(n!) - log(k!) - log(n-k!).
        /// </summary>
        /// <param name="n">Value.</param>
        /// <param name="k">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 LogBinomial(Complex32 n, Complex32 k)
        {
            return (Complex32)(GammaLog((Complex)n + 1) - GammaLog((Complex)k + 1) - GammaLog((Complex)n - (Complex)k + 1));
        }
        #endregion

        #region Laplace functions
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <param name="inverse">Reverse function or not.</param>
        /// <returns>Value.</returns>
        public static float Erf(float x, bool inverse)
        {
            return (float)(inverse ? InverseErf((double)x) : ErfValue((double)x));
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <param name="inverse">Reverse function or not.</param>
        /// <returns>Value.</returns>
        public static Complex32 Erf(Complex32 x, bool inverse)
        {
            return (Complex32)(inverse ? InverseErf((Complex)x) : ErfValue((Complex)x));
        }
        /// <summary>
        /// Returns the value of the imaginary error function.
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <returns>Value.</returns>
        public static float Erfi(float x)
        {
            if (Math.Abs(x) > 27) return x < 0 ? float.NegativeInfinity : float.PositiveInfinity;
            return (float)(-Complex.ImaginaryOne * ErfValue(Complex.ImaginaryOne * (double)x)).Real;
        }
        /// <summary>
        /// Returns the value of the imaginary error function.
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <returns>Value.</returns>
        public static Complex32 Erfi(Complex32 x)
        {
            return (Complex32)(-Complex.ImaginaryOne * ErfValue(Complex.ImaginaryOne * (Complex)x));
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <returns>Value.</returns>
        public static float Erf(float x)
        {
            return Erf(x, false);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <returns>Value.</returns>
        public static Complex32 Erf(Complex32 x)
        {
            return Erf(x, false);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <param name="a">The lower boundary of the normalization.</param>
        /// <param name="b">The upper limit of the normalization.</param>
        /// <returns>Value.</returns>
        public static float Erf(float x, float a, float b)
        {
            return (float)ErfValue(((double)x - (double)a) / (double)b);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <param name="a">The lower boundary of the normalization.</param>
        /// <param name="b">The upper limit of the normalization.</param>
        /// <returns>Value.</returns>
        public static Complex32 Erf(Complex32 x, Complex32 a, Complex32 b)
        {
            return (Complex32)ErfValue(((Complex)x - (Complex)a) / (Complex)b);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <returns>Value.</returns>
        public static float Erfc(float x)
        {
            return (float)ErfcValue((double)x);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <returns>Value.</returns>
        public static Complex32 Erfc(Complex32 x)
        {
            return (Complex32)ErfcValue((Complex)x);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <param name="a">The lower boundary of the normalization.</param>
        /// <param name="b">The upper limit of the normalization.</param>
        /// <returns>Value.</returns>
        public static float Erfc(float x, float a, float b)
        {
            return (float)ErfcValue(((double)x - (double)a) / (double)b);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral.</param>
        /// <param name="a">The lower boundary of the normalization.</param>
        /// <param name="b">The upper limit of the normalization.</param>
        /// <returns>Value.</returns>
        public static Complex32 Erfc(Complex32 x, Complex32 a, Complex32 b)
        {
            return (Complex32)ErfcValue(((Complex)x - (Complex)a) / (Complex)b);
        }



        #endregion

        #region Dawson function
        /// <summary>
        /// Returns the value of the D- / D + Dawson function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="positive">D- or D+.</param>
        /// <returns>Value.</returns>
        public static float Dawson(float x, bool positive)
        {
            if (!positive && Math.Abs(x) > 27) return x < 0 ? float.NegativeInfinity : float.PositiveInfinity;
            return (float)(positive ? DawsonValue((Complex)x) : -Complex.ImaginaryOne * DawsonValue(Complex.ImaginaryOne * (double)x)).Real;
        }
        /// <summary>
        /// Returns the value of the D- / D + Dawson function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="positive">D- or D+.</param>
        /// <returns>Value.</returns>
        public static Complex32 Dawson(Complex32 x, bool positive)
        {
            return (Complex32)(positive ? DawsonValue((Complex)x) : -Complex.ImaginaryOne * DawsonValue(Complex.ImaginaryOne * (Complex)x));
        }
        #endregion

        #region Faddeeva function
        /// <summary>
        /// Returns the value of the Faddeeva function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Faddeeva(float x)
        {
            return (Complex32)FaddeevaValue((Complex)x);
        }
        /// <summary>
        /// Returns the value of the Faddeeva function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Faddeeva(Complex32 x)
        {
            return (Complex32)FaddeevaValue((Complex)x);
        }
        #endregion

        #region Q-function
        /// <summary>
        /// Returns the value of a Q function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="inverse">Inverse function or not.</param>
        /// <returns>Value.</returns>
        public static float Q(float x, bool inverse = false)
        {
            return (float)(inverse ? Math.Sqrt(2) * InverseErfc(2.0 * x) : 0.5 * ErfcValue((double)x / Math.Sqrt(2)));
        }
        /// <summary>
        /// Returns the value of a Q function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="inverse">Inverse function or not.</param>
        /// <returns>Value.</returns>
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
        /// <see href="https://en.wikipedia.org/wiki/Hypergeometric_function"/>.
        /// </remarks>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <param name="c">Value.</param>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
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
        /// <see href="https://en.wikipedia.org/wiki/Hypergeometric_function"/>.
        /// </remarks>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <param name="c">Value.</param>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
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
        /// <see href="https://www.mathworks.com/help/symbolic/hypergeom.html#bt1nkmw-2"/>.
        /// </remarks>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
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
        /// <see href="https://www.mathworks.com/help/symbolic/hypergeom.html#bt1nkmw-2"/>.
        /// </remarks>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Hypergeom(Complex32 a, Complex32 b, Complex32 z)
        {
            return (Complex32)Hypergeometric1F1((Complex)a, (Complex)b, (Complex)z);
        }
        #endregion

        #region Beta functions
        /// <summary>
        /// Returns the value of the beta function: B(a, b) = Gamma(a) * Gamma(b) / Gamma(a + b).
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Value.</returns>
        public static float Beta(float a, float b)
        {
            return (float)BetaValue((double)a, (double)b);
        }
        /// <summary>
        /// Returns the value of the beta function: B(a, b) = Gamma(a) * Gamma(b) / Gamma(a + b).
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Beta(Complex32 a, Complex32 b)
        {
            return (Complex32)Complex.Exp(BetaLog((Complex)a, (Complex)b));
        }
        /// <summary>
        /// Returns the value of the beta function: B(m, n) = (m - 1)! * (n - 1)! / (m + n - 1)!.
        /// </summary>
        /// <param name="m">Integer number.</param>
        /// <param name="n">Integer number.</param>
        /// <returns>Value.</returns>
        public static double Beta(int m, int n)
        {
            return BetaValue((double)m, (double)n);
        }
        /// <summary>
        /// Returns the value of a derivative beta function: B'(a, b).
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Value.</returns>
        public static float BetaDerivative(float a, float b)
        {
            return (float)(Complex.Exp(BetaLog((Complex)a, (Complex)b)) * (Polygamma((Complex)a, false) - Polygamma((Complex)a + (Complex)b, false))).Real;
        }
        /// <summary>
        /// Returns the value of a derivative beta function: B'(a, b).
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 BetaDerivative(Complex32 a, Complex32 b)
        {
            return (Complex32)(Complex.Exp(BetaLog((Complex)a, (Complex)b)) * (Polygamma((Complex)a, false) - Polygamma((Complex)a + (Complex)b, false)));
        }
        /// <summary>
        /// Returns the value of an incomplete beta function: Bx(a, b).
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float BetaIncomplete(float a, float b, float x)
        {
            return (float)IncompleteBeta(a, b, x, false);
        }
        /// <summary>
        /// Returns the value of an incomplete beta function: Bx(a, b).
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
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
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float BetaIncompleteRegularized(float a, float b, float x)
        {
            return (float)IncompleteBeta(a, b, x, true);
        }
        /// <summary>
        /// Returns the value of a log-beta function.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Value.</returns>
        public static float LogBeta(float a, float b)
        {
            return (float)BetaLog((double)a, (double)b);
        }
        /// <summary>
        /// Returns the value of a log-beta function.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 LogBeta(Complex32 a, Complex32 b)
        {
            return (Complex32)BetaLog((Complex)a, (Complex)b);
        }
        #endregion

        #region Integral functions
        /// <summary>
        /// Returns the value of the integral cosine.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Ci(float x)
        {
            if (x < 0) return float.NaN;
            if (float.IsPositiveInfinity(x)) return 0;
            return (float)TrigonometricIntegral((Complex)x, false).Real;
        }
        /// <summary>
        /// Returns the value of the integral cosine.
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Ci(Complex32 z)
        {
            return (Complex32)TrigonometricIntegral((Complex)z, false);
        }
        /// <summary>
        /// Returns the value of the integral sine.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Si(float x)
        {
            if (float.IsInfinity(x)) return (float)(Math.Sign(x) * Math.PI / 2);
            return (float)TrigonometricIntegral((Complex)x, true).Real;
        }
        /// <summary>
        /// Returns the value of the integral sine.
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Si(Complex32 z)
        {
            return (Complex32)TrigonometricIntegral((Complex)z, true);
        }
        /// <summary>
        /// Returns the value of an integral exponential function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Ei(float x)
        {
            if (float.IsNegativeInfinity(x)) return 0;
            if (float.IsPositiveInfinity(x)) return float.PositiveInfinity;
            return (float)ExponentialIntegral((Complex)x).Real;
        }
        /// <summary>
        /// Returns the value of an integral exponential function.
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Ei(Complex32 z)
        {
            return (Complex32)ExponentialIntegral((Complex)z);
        }
        /// <summary>
        /// Returns the value of the integral logarithm.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Li(float x)
        {
            if (x < 0) return float.NaN;
            if (x == 0) return 0;
            return (float)ExponentialIntegral(Complex.Log((Complex)x)).Real;
        }
        /// <summary>
        /// Returns the value of the integral logarithm.
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
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
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float J(float x, int a)
        {
            return (float)BesselJ((Complex)x, a).Real;
        }
        /// <summary>
        /// Returns the value of a Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 J(Complex32 x, int a)
        {
            return (Complex32)BesselJ((Complex)x, a);
        }

        /// <summary>
        /// Returns the value of a Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Y(float x, int a)
        {
            if (x < 0) return float.NaN;
            return (float)BesselY((Complex)x, a).Real;
        }
        /// <summary>
        /// Returns the value of a Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Y(Complex32 x, int a)
        {
            return (Complex32)BesselY((Complex)x, a);
        }

        /// <summary>
        /// Returns the value of the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float I(float x, int a)
        {
            return (float)BesselI((Complex)x, a).Real;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 I(Complex32 x, int a)
        {
            return (Complex32)BesselI((Complex)x, a);
        }

        /// <summary>
        /// Returns the value of the modified Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float K(float x, int a)
        {
            if (x < 0) return float.NaN;
            return (float)BesselK((Complex)x, a).Real;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 K(Complex32 x, int a)
        {
            return (Complex32)BesselK((Complex)x, a);
        }



        #endregion

        #region Struve functions
        /// <summary>
        /// Returns the value of the Struve function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float H(float x, int a)
        {
            return (float)StruveValue((Complex)x, a, false).Real;
        }
        /// <summary>
        /// Returns the value of the Struve function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 H(Complex32 x, int a)
        {
            return (Complex32)StruveValue((Complex)x, a, false);
        }
        /// <summary>
        /// Returns the value of the modified Struve function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="v">Value.</param>
        /// <returns>Value.</returns>
        public static float L(float x, int v)
        {
            return (float)StruveValue((Complex)x, v, true).Real;
        }
        /// <summary>
        /// Returns the value of the modified Struve function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="v">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 L(Complex32 x, int v)
        {
            return (Complex32)StruveValue((Complex)x, v, true);
        }
        #endregion

        #region Fibonacci & Lucas numbers
        /// <summary>
        /// Returns the value of the Fibonacci number.
        /// </summary>
        /// <param name="n">Integer number.</param>
        /// <returns>Integer number.</returns>
        public static int Fibonacci(int n)
        {
            return FibonacciValue(n, false);
        }
        /// <summary>
        /// Returns the value of the Luca number.
        /// </summary>
        /// <param name="n">Integer number.</param>
        /// <returns>Integer number.</returns>
        public static int Lucas(int n)
        {
            return FibonacciValue(n, true);
        }
        #endregion

        #region Harmonic number
        /// <summary>
        /// Returns the harmonic number.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <returns>Value.</returns>
        public static float Harm(int n)
        {
            if (n < 0) return float.NaN;
            if (n == 0) return 0;
            return (float)(Polygamma((Complex)((double)n + 1), false).Real + EulerGamma);
        }
        /// <summary>
        /// Returns the harmonic number.
        /// </summary>
        /// <param name="n">Order.</param>
        /// <param name="m">Value.</param>
        /// <returns>Value.</returns>
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
        /// <param name="n">Value.</param>
        /// <returns>Value.</returns>
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
        /// <param name="n">Order.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static double Euler(int n, float x)
        {
            return NumberPolynomial(n, x, true);
        }
        #endregion

        #region Bernoulli function
        /// <summary>
        /// Returns the Bernoulli number.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <returns>Value.</returns>
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
        /// <param name="n">Order.</param>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static double Bernoulli(int n, float x)
        {
            return NumberPolynomial(n, x, false);
        }
        #endregion

        #region Minkowski function
        /// <summary>
        /// Returns the value of the Minkowski function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
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

        #region Private gamma, beta, and recurrence kernels
        // Keep intermediate arithmetic in double precision. Public single-precision
        // overloads round only once, after the tail or logarithmic ratio is evaluated.
        /// <summary>
        /// Relative stopping threshold for double-precision series and continued fractions.
        /// </summary>
        private const double SpecialEpsilon = 2e-15;
        /// <summary>
        /// Square root of pi used by error-function and Bessel kernels.
        /// </summary>
        private const double SqrtPi = 1.7724538509055160273;
        /// <summary>
        /// Complex sentinel for undefined values or unsuccessful numerical convergence.
        /// </summary>
        private static readonly Complex ComplexNaN = new Complex(double.NaN, double.NaN);
        /// <summary>
        /// Lanczos coefficients for the gamma approximation with shift g = 7.
        /// </summary>
        private static readonly double[] GammaCoefficients =
        {
            0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
        };

        /// <summary>
        /// Tests whether the argument is a nonpositive integer pole of the gamma function.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>True exactly at a nonpositive real integer.</returns>
        private static bool IsPole(double x) => x <= 0 && x == Math.Floor(x);
        /// <summary>
        /// Tests whether the argument is a nonpositive integer pole of the gamma function.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>True exactly at a nonpositive real integer.</returns>
        private static bool IsPole(Complex x) => x.Imaginary == 0 && IsPole(x.Real);
        /// <summary>
        /// Tests whether both components of a complex argument are finite.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>True if neither component is NaN or infinite.</returns>
        private static bool IsFinite(Complex x) => !double.IsNaN(x.Real) && !double.IsInfinity(x.Real)
            && !double.IsNaN(x.Imaginary) && !double.IsInfinity(x.Imaginary);

        /// <summary>
        /// Evaluates log-gamma using Lanczos coefficients, recurrence, and reflection.
        /// </summary>
        /// <remarks>The real overload returns log(abs(Gamma(x))). The complex overload retains the analytic logarithm on the plane cut along the negative real axis.</remarks>
        /// <param name="x">Function argument.</param>
        /// <returns>The logarithm of the gamma magnitude for real inputs, or the analytic log-gamma for complex inputs.</returns>
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

        /// <summary>
        /// Returns the sign of gamma at a real argument away from its poles.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>One or minus one according to the sign of Gamma(x).</returns>
        private static double GammaSign(double x) => x > 0 ? 1 : Math.Sign(Math.Sin(Math.PI * (x % 2)));

        /// <summary>
        /// Evaluates gamma from its logarithmic magnitude or analytic logarithm.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>Gamma at the argument, or NaN at a pole or an undefined input.</returns>
        private static double GammaValue(double x)
        {
            double log = GammaLog(x);
            return double.IsNaN(log) ? double.NaN : GammaSign(x) * Math.Exp(log);
        }

        // Analytic log-gamma on the plane cut along the negative real axis.
        // Recurrence retains the winding of log Gamma instead of taking log(Gamma).
        /// <summary>
        /// Evaluates log-gamma using Lanczos coefficients, recurrence, and reflection.
        /// </summary>
        /// <remarks>The real overload returns log(abs(Gamma(x))). The complex overload retains the analytic logarithm on the plane cut along the negative real axis.</remarks>
        /// <param name="x">Function argument.</param>
        /// <returns>The logarithm of the gamma magnitude for real inputs, or the analytic log-gamma for complex inputs.</returns>
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

        /// <summary>
        /// Evaluates gamma from its logarithmic magnitude or analytic logarithm.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>Gamma at the argument, or NaN at a pole or an undefined input.</returns>
        private static Complex GammaValue(Complex x) => Complex.Exp(GammaLog(x));

        /// <summary>
        /// Evaluates the lower incomplete-gamma series before multiplication by x^s*exp(-x).
        /// </summary>
        /// <param name="s">Gamma shape parameter away from series denominator poles.</param>
        /// <param name="x">Integration argument in the region selected for series evaluation.</param>
        /// <returns>The lower-gamma series factor, or NaN if the iteration limit is reached.</returns>
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

        /// <summary>
        /// Evaluates the lower incomplete-gamma series before multiplication by x^s*exp(-x).
        /// </summary>
        /// <param name="s">Gamma shape parameter away from series denominator poles.</param>
        /// <param name="x">Integration argument in the region selected for series evaluation.</param>
        /// <returns>The lower-gamma series factor, or NaN if the iteration limit is reached.</returns>
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

        /// <summary>
        /// Evaluates the upper incomplete-gamma continued-fraction factor before multiplication by x^s*exp(-x).
        /// </summary>
        /// <param name="s">Gamma shape parameter.</param>
        /// <param name="x">Integration argument in the region selected for continued-fraction evaluation.</param>
        /// <returns>The upper-gamma continued-fraction factor, or NaN if the iteration limit is reached.</returns>
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

        /// <summary>
        /// Evaluates the upper incomplete-gamma continued-fraction factor before multiplication by x^s*exp(-x).
        /// </summary>
        /// <param name="s">Gamma shape parameter.</param>
        /// <param name="x">Integration argument in the region selected for continued-fraction evaluation.</param>
        /// <returns>The upper-gamma continued-fraction factor, or NaN if the iteration limit is reached.</returns>
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

        /// <summary>
        /// Evaluates a lower or upper incomplete gamma function, optionally divided by gamma(s).
        /// </summary>
        /// <param name="s">Shape parameter; real overloads require a positive finite value.</param>
        /// <param name="x">Integration limit; real overloads require a nonnegative value.</param>
        /// <param name="upper">True selects the upper tail; false selects the lower integral.</param>
        /// <param name="regularized">True divides by Gamma(s).</param>
        /// <returns>The requested incomplete gamma value or regularized ratio, or NaN for invalid inputs or nonconvergence.</returns>
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

        /// <summary>
        /// Evaluates a lower or upper incomplete gamma function, optionally divided by gamma(s).
        /// </summary>
        /// <param name="s">Shape parameter; real overloads require a positive finite value.</param>
        /// <param name="x">Integration limit; real overloads require a nonnegative value.</param>
        /// <param name="upper">True selects the upper tail; false selects the lower integral.</param>
        /// <param name="regularized">True divides by Gamma(s).</param>
        /// <returns>The requested incomplete gamma value or regularized ratio, or NaN for invalid inputs or nonconvergence.</returns>
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

        /// <summary>
        /// Evaluates the logarithmic beta ratio from three log-gamma values.
        /// </summary>
        /// <param name="a">First beta parameter.</param>
        /// <param name="b">Second beta parameter.</param>
        /// <returns>The log-gamma combination logGamma(a) + logGamma(b) - logGamma(a + b).</returns>
        private static double BetaLog(double a, double b) => GammaLog(a) + GammaLog(b) - GammaLog(a + b);
        /// <summary>
        /// Evaluates the logarithmic beta ratio from three log-gamma values.
        /// </summary>
        /// <param name="a">First beta parameter.</param>
        /// <param name="b">Second beta parameter.</param>
        /// <returns>The log-gamma combination logGamma(a) + logGamma(b) - logGamma(a + b).</returns>
        private static Complex BetaLog(Complex a, Complex b) => GammaLog(a) + GammaLog(b) - GammaLog(a + b);
        /// <summary>
        /// Evaluates real beta with gamma-sign tracking and a logarithmic magnitude.
        /// </summary>
        /// <param name="a">First real beta parameter.</param>
        /// <param name="b">Second real beta parameter.</param>
        /// <returns>The signed real beta value, or NaN where the gamma ratio is undefined.</returns>
        private static double BetaValue(double a, double b)
        {
            double log = BetaLog(a, b);
            return double.IsNaN(log) ? double.NaN : GammaSign(a) * GammaSign(b) * GammaSign(a + b) * Math.Exp(log);
        }

        /// <summary>
        /// Evaluates the continued-fraction factor used by the incomplete beta function.
        /// </summary>
        /// <param name="a">Positive first beta shape parameter.</param>
        /// <param name="b">Positive second beta shape parameter.</param>
        /// <param name="x">Integration limit in (0, 1).</param>
        /// <returns>The incomplete-beta continued-fraction factor, or NaN if the iteration limit is reached.</returns>
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

        /// <summary>
        /// Evaluates the lower incomplete beta function, using a complementary tail when it improves conditioning.
        /// </summary>
        /// <param name="a">Positive finite first beta shape parameter.</param>
        /// <param name="b">Positive finite second beta shape parameter.</param>
        /// <param name="x">Integration limit in [0, 1].</param>
        /// <param name="regularized">True divides by B(a, b).</param>
        /// <returns>The lower incomplete beta value or regularized ratio, or NaN outside the supported domain.</returns>
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

        /// <summary>
        /// Evaluates a first- or second-kind Chebyshev polynomial, including exact endpoint values and negative-order identities.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <param name="n">Integer order.</param>
        /// <param name="secondKind">True selects the second kind; false selects the first kind.</param>
        /// <returns>T_n(x) or U_n(x); Int32.MinValue returns NaN.</returns>
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

        /// <summary>
        /// Evaluates a rising or falling factorial using finite products or a continued gamma ratio.
        /// </summary>
        /// <param name="n">Complex initial factor.</param>
        /// <param name="k">Complex product order; nonnegative integer orders use a finite product.</param>
        /// <param name="rising">True selects the rising factorial; false selects the falling factorial.</param>
        /// <returns>The rising or falling factorial, continued through the log-gamma ratio for noninteger orders.</returns>
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

        /// <summary>
        /// Evaluates a generalized binomial coefficient using finite products or logarithmic gamma ratios.
        /// </summary>
        /// <param name="n">Complex upper index.</param>
        /// <param name="k">Complex lower index.</param>
        /// <returns>The generalized binomial coefficient; a negative real lower index returns zero.</returns>
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

        /// <summary>
        /// Evaluates Gauss 2F1 using terminating series, transformations, or differential-equation continuation.
        /// </summary>
        /// <param name="a">First numerator parameter.</param>
        /// <param name="b">Second numerator parameter.</param>
        /// <param name="c">Denominator parameter; nonpositive real integers are rejected.</param>
        /// <param name="z">Complex function argument.</param>
        /// <returns>The continued 2F1 value, or NaN at a rejected parameter pole or failed convergence.</returns>
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
        #endregion

        #region Private elementary and number-polynomial kernels
        /// <summary>
        /// Evaluates exp(x) - 1 with a power series near zero to preserve small results.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>exp(x) - 1, retaining small differences near zero.</returns>
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

        /// <summary>
        /// Evaluates 1 / (1 + exp(-z)), selecting the exponential with a nonpositive real part.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <returns>The complex logistic value 1 / (1 + exp(-z)).</returns>
        private static Complex LogisticValue(Complex z)
        {
            if (z.Real >= 0) return 1 / (1 + Complex.Exp(-z));
            Complex exponential = Complex.Exp(z);
            return exponential / (1 + exponential);
        }

        /// <summary>
        /// Evaluates the normalized sign of sin(2^n*pi*z), using zero at real integer nodes.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="n">Integer order.</param>
        /// <returns>The normalized complex sine, with zero at real nodes and NaN for nonfinite inputs.</returns>
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

        /// <summary>
        /// Computes a Fibonacci or Lucas number exactly over the supported signed Int32 result range.
        /// </summary>
        /// <param name="n">Signed sequence index: [-46, 46] for Fibonacci or [-44, 44] for Lucas.</param>
        /// <param name="lucas">True selects Lucas numbers; false selects Fibonacci numbers.</param>
        /// <returns>The exact signed sequence value; out-of-range requests throw ArgumentOutOfRangeException.</returns>
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

        /// <summary>
        /// Evaluates a Bernoulli or Euler polynomial by Horner accumulation of number-table coefficients.
        /// </summary>
        /// <param name="n">Polynomial degree: 0 through 186 for Euler or 258 for Bernoulli.</param>
        /// <param name="x">Function argument.</param>
        /// <param name="euler">True selects Euler polynomials; false selects Bernoulli polynomials.</param>
        /// <returns>The selected polynomial value, or NaN for an unsupported order.</returns>
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
        #endregion

        #region Private series and analytic-continuation kernels
        /// <summary>
        /// Evaluates a classical orthogonal polynomial by its three-term recurrence.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <param name="a">Laguerre alpha or Gegenbauer lambda parameter; unused for Legendre and Hermite.</param>
        /// <param name="n">Nonnegative polynomial degree.</param>
        /// <param name="family">Polynomial family: 0 Laguerre, 1 Legendre, 2 Hermite, otherwise Gegenbauer.</param>
        /// <returns>The selected degree-n polynomial at x, or NaN for negative degree.</returns>
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

        /// <summary>
        /// Evaluates Erlang B blocking probability by recurrence, avoiding direct powers and factorials.
        /// </summary>
        /// <param name="y">Offered traffic argument.</param>
        /// <param name="n">Nonnegative number of servers.</param>
        /// <returns>The Erlang B value, or NaN for negative server counts or nonfinite traffic.</returns>
        private static Complex ErlangBlocking(Complex y, int n)
        {
            if (n < 0 || !IsFinite(y)) return ComplexNaN;
            Complex blocking = 1;
            for (int k = 1; k <= n; k++) blocking = y * blocking / (k + y * blocking);
            return blocking;
        }

        /// <summary>
        /// Evaluates the principal log(sin(pi*z)) without overflowing at large imaginary arguments.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <returns>The natural logarithm of the sine magnitude and its principal phase.</returns>
        private static Complex LogSinPi(Complex z)
        {
            double real = z.Real % 2;
            if (Math.Abs(z.Imaginary) < 20) return Complex.Log(Complex.Sin(Math.PI * new Complex(real, z.Imaginary)));
            double phase = (z.Imaginary > 0 ? 1 : -1) * (Math.PI / 2 - Math.PI * real);
            return new Complex(Math.PI * Math.Abs(z.Imaginary) - Math.Log(2), Math.Atan2(Math.Sin(phase), Math.Cos(phase)));
        }

        /// <summary>
        /// Evaluates digamma or trigamma by reflection, recurrence, and a Bernoulli asymptotic expansion.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="derivative">True selects trigamma; false selects digamma.</param>
        /// <returns>Digamma or trigamma at z, or NaN at poles or nonfinite arguments.</returns>
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

        /// <summary>
        /// Evaluates the Riemann zeta function using reflection and Euler–Maclaurin summation.
        /// </summary>
        /// <param name="s">Complex zeta argument; s = 1 is a pole.</param>
        /// <returns>Zeta(s), or NaN at the pole s = 1 or for nonfinite inputs.</returns>
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
        /// <summary>
        /// Advances a hypergeometric function and its first derivative by a Taylor step away from singularities.
        /// </summary>
        /// <param name="a">First numerator parameter; unused for 0F1.</param>
        /// <param name="b">Second numerator parameter for 2F1, denominator parameter otherwise.</param>
        /// <param name="c">Denominator parameter for 2F1; unused otherwise.</param>
        /// <param name="center">Current center of the Taylor expansion.</param>
        /// <param name="step">Nonzero complex step within the convergence disk.</param>
        /// <param name="kind">Hypergeometric family: 0 for 0F1, 1 for 1F1, 2 for 2F1.</param>
        /// <param name="value">On entry, the value at center; on return, the value at center + step.</param>
        /// <param name="derivative">On entry, the first derivative at center; on return, the derivative at center + step.</param>
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

        /// <summary>
        /// Continues a hypergeometric solution through overlapping Taylor disks that avoid the differential equation's singularities.
        /// </summary>
        /// <param name="a">First numerator parameter; unused for 0F1.</param>
        /// <param name="b">Second numerator parameter for 2F1, denominator parameter otherwise.</param>
        /// <param name="c">Denominator parameter for 2F1; unused otherwise.</param>
        /// <param name="target">Target point for analytic continuation.</param>
        /// <param name="kind">Hypergeometric family: 0 for 0F1, 1 for 1F1, 2 for 2F1.</param>
        /// <returns>The solution at target, or NaN if the path cannot be completed within the iteration limit.</returns>
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

        /// <summary>
        /// Evaluates confluent hypergeometric functions, including reduced families selected by absent-parameter sentinels.
        /// </summary>
        /// <param name="a">Numerator parameter; NaN denotes an absent numerator parameter.</param>
        /// <param name="b">Denominator parameter; NaN denotes an absent denominator parameter.</param>
        /// <param name="z">Complex function argument.</param>
        /// <returns>The 1F1, 0F1, 1F0, or 0F0 value selected by the supplied parameters.</returns>
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
        #endregion

        #region Private error-function and quadrature kernels
        /// <summary>
        /// Positive nodes of the symmetric 16-point Gauss–Legendre rule on [-1, 1].
        /// </summary>
        private static readonly double[] GaussNodes =
        {
            0.0950125098376374402, 0.281603550779258913, 0.458016777657227386, 0.617876244402643748,
            0.755404408355003034, 0.865631202387831744, 0.944575023073232576, 0.989400934991649933
        };
        /// <summary>
        /// Weights paired with the positive and negative Gauss–Legendre nodes.
        /// </summary>
        private static readonly double[] GaussWeights =
        {
            0.189450610455068496, 0.182603415044923589, 0.169156519395002538, 0.149595988816576733,
            0.124628971255533872, 0.0951585116824927848, 0.0622535239386478929, 0.0271524594117540949
        };

        /// <summary>
        /// Integrates a complex-valued function over a real interval using composite 16-point Gauss–Legendre quadrature.
        /// </summary>
        /// <param name="f">Integrand evaluated at real quadrature nodes.</param>
        /// <param name="end">Real upper integration limit; the lower limit is zero.</param>
        /// <param name="panels">Positive number of equal quadrature panels.</param>
        /// <returns>The quadrature estimate of the integral from zero to end.</returns>
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
        /// <summary>
        /// Evaluates exp(-z^2)*erfc(-i*z) using reflection, decaying quadrature, and an asymptotic expansion.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <returns>The Faddeeva function w(z), or NaN for nonfinite inputs.</returns>
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

        /// <summary>
        /// Evaluates the complementary error function directly so small tails are not obtained by subtracting from one.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>The complementary error function at the argument.</returns>
        private static double ErfcValue(double x)
        {
            if (double.IsNaN(x)) return double.NaN;
            if (x < 0) return 2 - ErfcValue(-x);
            if (double.IsPositiveInfinity(x)) return 0;
            return IncompleteGamma(0.5, x * x, true, true);
        }

        /// <summary>
        /// Evaluates the error function using its local power series and complementary-function relations.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>The error function at the argument.</returns>
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

        /// <summary>
        /// Evaluates the complementary error function directly so small tails are not obtained by subtracting from one.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <returns>The complementary error function at the argument.</returns>
        private static Complex ErfcValue(Complex z)
        {
            if (z.Imaginary == 0) return new Complex(ErfcValue(z.Real), 0);
            if (z.Real < 0) return 2 - ErfcValue(-z);
            return Complex.Exp(-z * z) * FaddeevaValue(Complex.ImaginaryOne * z);
        }

        /// <summary>
        /// Evaluates the error function using its local power series and complementary-function relations.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <returns>The error function at the argument.</returns>
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

        /// <summary>
        /// Inverts the real complementary error function by bracketing, including explicit endpoint limits.
        /// </summary>
        /// <param name="p">Complementary-error-function value in [0, 2].</param>
        /// <returns>The inverse erfc, with +infinity at zero, -infinity at two, and NaN outside [0, 2].</returns>
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

        /// <summary>
        /// Inverts the error function, using real-domain limits or complex Halley refinement as appropriate.
        /// </summary>
        /// <param name="x">Real error-function value in [-1, 1].</param>
        /// <returns>An inverse-error-function value; invalid real inputs or failed complex refinement return NaN.</returns>
        private static double InverseErf(double x)
        {
            if (double.IsNaN(x) || Math.Abs(x) > 1) return double.NaN;
            if (Math.Abs(x) < 1e-4) return SqrtPi / 2 * x * (1 + Math.PI * x * x / 12);
            return x < 0 ? -InverseErfc(1 + x) : InverseErfc(1 - x);
        }

        /// <summary>
        /// Inverts the error function, using real-domain limits or complex Halley refinement as appropriate.
        /// </summary>
        /// <param name="z">Complex error-function value to invert.</param>
        /// <returns>An inverse-error-function value; invalid real inputs or failed complex refinement return NaN.</returns>
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

        /// <summary>
        /// Evaluates the entire continuation of Dawson's integral using a local series and the Faddeeva function.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <returns>exp(-z^2) times the integral of exp(t^2) from zero to z.</returns>
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

        /// <summary>
        /// Evaluates the library's unnormalized Fresnel integral of sin(t^2) or cos(t^2) from zero to z.
        /// </summary>
        /// <remarks>The integrand uses t^2, without the pi/2 factor of the normalized Fresnel convention.</remarks>
        /// <param name="z">Complex function argument.</param>
        /// <param name="sine">True selects the sine integral; false selects the cosine integral.</param>
        /// <returns>The integral of the selected unnormalized Fresnel integrand from zero to z.</returns>
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

        /// <summary>
        /// Evaluates the entire continuation of n!/sqrt(pi) times the integral of exp(-t^n) from zero to z.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="n">Nonnegative integer exponent in exp(-t^n).</param>
        /// <returns>The generalized error integral, or NaN for negative order or nonfinite z.</returns>
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

        /// <summary>
        /// Solves w*exp(w) = z on the requested Lambert W branch using a branch-aware initial value and refinement.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="branch">Integer Lambert W branch index.</param>
        /// <returns>W_branch(z), or NaN for rejected inputs or failed refinement.</returns>
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
        #endregion

        #region Private special-integral kernels
        /// <summary>
        /// Evaluates E1(z) on its logarithmic branch using a local series or an upper-gamma continued fraction.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <returns>The continued E1 value using the principal logarithm in its local series.</returns>
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

        /// <summary>
        /// Evaluates Ei(z) with explicit continuation terms around the negative real axis.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <returns>Ei(z), real on the negative real axis and with signed imaginary continuation off that axis.</returns>
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

        /// <summary>
        /// Evaluates Si(z) or Ci(z) using local series and exponential-integral relations.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="sine">True selects the sine integral; false selects the cosine integral.</param>
        /// <returns>Si(z) or Ci(z); the cosine integral uses the upper value on the negative real axis.</returns>
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

        /// <summary>
        /// Evaluates Owen's T function by quadrature after the substitution t = tan(theta).
        /// </summary>
        /// <param name="h">Gaussian threshold, continued to complex values.</param>
        /// <param name="a">Integration limit, continued through the principal complex arctangent.</param>
        /// <returns>The continued Owen T value from the transformed quadrature.</returns>
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
        #endregion

        #region Private Bessel and Struve kernels
        /// <summary>
        /// Euler–Mascheroni constant.
        /// </summary>
        private const double EulerGamma = 0.57721566490153286061;

        /// <summary>
        /// Returns i raised to an integer power by reducing the exponent modulo four.
        /// </summary>
        /// <param name="n">Signed integer exponent.</param>
        /// <returns>One of 1, i, -1, or -i.</returns>
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

        /// <summary>
        /// Evaluates the power series for Bessel J or modified Bessel I at a nonnegative integer order.
        /// </summary>
        /// <param name="z">Nonzero complex argument in the series region.</param>
        /// <param name="n">Nonnegative integer order.</param>
        /// <param name="modified">True selects the modified function; false selects the ordinary function.</param>
        /// <returns>J_n(z) or I_n(z), or NaN if the series iteration limit is reached.</returns>
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
        /// <summary>
        /// Evaluates the Hankel asymptotic expansion for Bessel J or Y of order zero or one.
        /// </summary>
        /// <remarks>Used only for base orders zero and one; higher orders use recurrence.</remarks>
        /// <param name="z">Nonzero argument in the asymptotic region.</param>
        /// <param name="n">Base order, zero or one.</param>
        /// <param name="secondKind">True selects the second kind; false selects the first kind.</param>
        /// <returns>The asymptotic approximation to J_n(z) or Y_n(z).</returns>
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

        /// <summary>
        /// Evaluates integer-order Bessel J using a power series, asymptotics, and stable upward or Miller recurrence.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="order">Integer order in [-100000, 100000].</param>
        /// <returns>J_order(z), or NaN for a rejected order or nonfinite input.</returns>
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

        /// <summary>
        /// Evaluates Bessel Y of order zero or one from its logarithmic series or asymptotic expansion.
        /// </summary>
        /// <param name="z">Nonzero argument; the caller uses positive real inputs.</param>
        /// <param name="n">Base order, zero or one.</param>
        /// <returns>Y_n(z) for the selected base order.</returns>
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

        /// <summary>
        /// Evaluates integer-order Bessel Y, retaining the upper/lower continuation across the negative real axis.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="order">Integer order in [-100000, 100000].</param>
        /// <returns>Y_order(z), with the upper continuation on the negative real axis.</returns>
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

        /// <summary>
        /// Evaluates integer-order modified Bessel I by rotating the argument of Bessel J.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="n">Integer order in [-100000, 100000].</param>
        /// <returns>I_n(z), using I_(-n)(z) = I_n(z).</returns>
        private static Complex BesselI(Complex z, int n)
        {
            if (n == int.MinValue) return ComplexNaN;
            n = Math.Abs(n);
            return ImaginaryPower(-n) * BesselJ(Complex.ImaginaryOne * z, n);
        }

        /// <summary>
        /// Evaluates modified Bessel K of order zero or one using a logarithmic series or a decaying asymptotic expansion.
        /// </summary>
        /// <param name="z">Nonzero complex argument in the closed right half-plane.</param>
        /// <param name="n">Base order, zero or one.</param>
        /// <returns>K_n(z) for the selected base order.</returns>
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

        /// <summary>
        /// Evaluates integer-order modified Bessel K using base orders, recurrence, and continuation across its cut.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="order">Integer order in [-100000, 100000].</param>
        /// <returns>K_order(z), with the upper continuation on the negative real axis.</returns>
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

        /// <summary>
        /// Evaluates ordinary Struve H or modified Struve L using a series, recurrence, or transformed quadrature.
        /// </summary>
        /// <param name="z">Complex function argument.</param>
        /// <param name="n">Integer order in [-100000, 100000].</param>
        /// <param name="modified">True selects the modified function; false selects the ordinary function.</param>
        /// <returns>H_n(z) or L_n(z), or NaN for rejected inputs.</returns>
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
        #endregion

        #region Internal probability-distribution kernels
        // Reuse the audited double kernels without rounding intermediate distribution values to float.
        /// <summary>
        /// Provides the double-precision log-gamma kernel for probability-distribution calculations.
        /// </summary>
        /// <param name="x">Positive real gamma argument.</param>
        /// <returns>The natural logarithm of Gamma(x) for positive x.</returns>
        internal static double DistributionLogGamma(double x) => GammaLog(x);
        /// <summary>
        /// Provides the double-precision log-beta kernel for positive distribution shape parameters.
        /// </summary>
        /// <param name="a">Positive first beta shape parameter.</param>
        /// <param name="b">Positive second beta shape parameter.</param>
        /// <returns>The natural logarithm of B(a, b).</returns>
        internal static double DistributionLogBeta(double a, double b) => BetaLog(a, b);
        /// <summary>
        /// Evaluates the real digamma function for positive distribution shape parameters.
        /// </summary>
        /// <param name="x">Positive real gamma argument.</param>
        /// <returns>The logarithmic derivative of Gamma(x).</returns>
        internal static double DistributionDigamma(double x) => Polygamma(new Complex(x, 0), false).Real;
        /// <summary>
        /// Evaluates a regularized lower or upper incomplete gamma function without intermediate float rounding.
        /// </summary>
        /// <param name="a">Positive finite gamma shape parameter.</param>
        /// <param name="x">Nonnegative integration limit.</param>
        /// <param name="upper">True selects the upper tail; false selects the lower integral.</param>
        /// <returns>P(a, x) or Q(a, x), according to upper.</returns>
        internal static double DistributionGamma(double a, double x, bool upper) => IncompleteGamma(a, x, upper, true);
        /// <summary>
        /// Evaluates the regularized incomplete beta function without intermediate float rounding.
        /// </summary>
        /// <param name="a">Positive finite first beta shape parameter.</param>
        /// <param name="b">Positive finite second beta shape parameter.</param>
        /// <param name="x">Integration limit in [0, 1].</param>
        /// <returns>The regularized lower ratio I_x(a, b).</returns>
        internal static double DistributionBeta(double a, double b, double x) => IncompleteBeta(a, b, x, true);
        /// <summary>
        /// Evaluates the real complementary error function in double precision for probability tails.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>erfc(x), including the real infinite-argument limits.</returns>
        internal static double DistributionErfc(double x) => ErfcValue(x);
        /// <summary>
        /// Evaluates exp(x) - 1 without cancellation near zero for distribution calculations.
        /// </summary>
        /// <param name="x">Function argument.</param>
        /// <returns>exp(x) - 1, retaining small differences near zero.</returns>
        internal static double DistributionExpm1(double x) => Expm1(x);
        #endregion
    }
}
