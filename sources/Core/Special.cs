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
    public static class Special
    {
        #region Private data
        /// <summary>
        /// 1 / ( 2 * PI ).
        /// </summary>
        private const float INV_2PI = 0.15915494309189533576888376337251f;
        /// <summary>
        /// Sqrt( PI ).
        /// </summary>
        private const float SQRT_PI = 1.7724538509055160272981674833411f;
        #endregion

        #region Chebyshev polynomial
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static float ChebyshevT(float x, int n)
        {
            return Maths.Cos(n * Maths.Acos(x));
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static Complex32 ChebyshevT(Complex32 x, int n)
        {
            return Maths.Cos(n * Maths.Acos(x));
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static float ChebyshevU(float x, int n)
        {
            float z = Maths.Acos(x);
            return Maths.Sin((n + 1) * z) / Maths.Sin(z);
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static Complex32 ChebyshevU(Complex32 x, int n)
        {
            Complex32 z = Maths.Acos(x);
            return Maths.Sin((n + 1) * z) / Maths.Sin(z);
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
            if (n == 1) return x;

            // Generalized formula
            // Abel polynomials recurrence relation for any n ≥ 1:
            return x * Maths.Pow(x - a * n, n - 1);
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
            if (n == 1) return x;

            // Generalized formula
            // Abel polynomials recurrence relation for any n ≥ 1:
            return x * Maths.Pow(x - a * n, n - 1);
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
            if (k < 0) return float.NaN;
            if (k == 0) return 1f;
            if (k == 1) return 1f + a - x;

            float num = (2f * k - 1f + a - x) * Laguerre(x, a, k - 1) - (k - 1f + a) * Laguerre(x, a, k - 2);
            return num / k;
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
            if (k < 0) return Complex32.NaN;
            if (k == 0) return 1f;
            if (k == 1) return 1f + a - x;

            Complex32 num = (2f * k - 1f + a - x) * Laguerre(x, a, k - 1) - (k - 1f + a) * Laguerre(x, a, k - 2);
            return num / k;
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
            if (m < 0) return float.NaN;
            if (m == 0) return 1f;
            if (m == 1) return x;

            float num = (2f * m - 1f) * x * Legendre(x, m - 1) - (m - 1f) * Legendre(x, m - 2);
            return num / m;
        }
        /// <summary>
        /// Returns the value of the Legendre polynomial of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="m">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Legendre(Complex32 x, int m)
        {
            if (m < 0) return Complex32.NaN;
            if (m == 0) return 1f;
            if (m == 1) return x;

            Complex32 num = (2f * m - 1f) * x * Legendre(x, m - 1) - (m - 1f) * Legendre(x, m - 2);
            return num / m;
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
            if (m < 0) return float.NaN;
            if (m == 0) return 1f;
            if (m == 1) return 2f * x;

            float ksi = x * Hermite(x, m - 1) - (m - 1) * Hermite(x, m - 2);
            return 2f * ksi;
        }
        /// <summary>
        /// Returns the value of the Hermite polynomial.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="m">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Hermite(Complex32 x, int m)
        {
            if (m < 0) return Complex32.NaN;
            if (m == 0) return Complex32.One;
            if (m == 1) return 2f * x;

            Complex32 ksi = x * Hermite(x, m - 1) - (m - 1) * Hermite(x, m - 2);
            return 2f * ksi;
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
            if (n < 0) return float.NaN;
            if (n == 0) return 1;
            if (n == 1) return 2 * a * x;

            // Generalized formula
            // Laguerre polynomials recurrence relation for any k ≥ 1:
            float psi = 2.0f * x * (n + a - 1) * Gegenbauer(x, a, n - 1) - (n + 2 * a - 2) * Gegenbauer(x, a, n - 2);
            return psi / n;
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
            if (n < 0) return Complex32.NaN;
            if (n == 0) return 1;
            if (n == 1) return 2 * a * x;

            // Generalized formula
            // Laguerre polynomials recurrence relation for any k ≥ 1:
            Complex32 psi = 2.0f * x * (n + a - 1) * Gegenbauer(x, a, n - 1) - (n + 2 * a - 2) * Gegenbauer(x, a, n - 2);
            return psi / n;
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
            return Maths.Asin(Maths.Tanh(x));
        }
        /// <summary>
        /// Returns the value of the Guderman function.
        /// </summary>
        /// <param name="x">Angle in radians</param>
        /// <returns>Value</returns>
        public static Complex32 Gd(Complex32 x)
        {
            return Maths.Asin(Maths.Tanh(x));
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
            float p = Maths.Pow(2, n);
            float v = p * Maths.Pi * t;
            return Math.Sign(Math.Sin(v));
        }
        /// <summary>
        /// Returns the value of the Radamecher function.
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="n">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Rademacher(Complex32 z, int n)
        {
            float p = Maths.Pow(2.0f, n);
            Complex32 v = new Complex32(p * Maths.Pi, 0f) * z;
            Complex32 s = Maths.Sin(v);

            float mag = Maths.Abs(s);
            if (mag == 0f) return Complex32.Zero;  // zeros at z = k / 2^n for real z

            return s / mag; // complex signum
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
            return 0.5f + 0.5f * Maths.Tanh(k * x);
        }
        /// <summary>
        /// Returns the value of the Heaviside delta function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="k">Smoothing factor</param>
        /// <returns>Value</returns>
        public static Complex32 Heaviside(Complex32 x, Complex32 k)
        {
            return 0.5f + 0.5f * Maths.Tanh(k * x);
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
            float x = -c * t;
            float y = -b * Maths.E;
            float z = a * Maths.E;
            return Maths.Pow(z, Maths.Pow(y, x));
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
            Complex32 x = -c * t;
            Complex32 y = -b * Maths.E;
            Complex32 z = a * Maths.E;
            return Maths.Pow(z, Maths.Pow(y, x));
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
            return Special.Logistic(x, a, k, b, 1.0f, 1.0f, 1.0f);
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
            return Special.Logistic(x, a, k, b, 1.0f, 1.0f, 1.0f);
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Logistic(float x)
        {
            return Special.Logistic(x, 0, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Logistic(Complex32 x)
        {
            return Special.Logistic(x, 0, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
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
            float e = Special.Erlang(y, v);
            float a = v * e;
            float b = v - y + y * e;
            float c = (v - y) * t;
            return a / b * Maths.Exp(-c);
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
            Complex32 e = Special.Erlang(y, v);             // Erlang-B (blocking)
            Complex32 a = v * e;                    // v * B
            Complex32 b = v - y + y * e;            // v - y + y B = v(1-ρ) + ρ v B
            Complex32 c = (v - y) * t;              // (v - y) t
            return a / b * Maths.Exp(-c);           // C * exp( - (v - y) t )
        }
        /// <summary>
        /// Returns the value of the Erlang B-function.
        /// </summary>
        /// <param name="y">First parameter</param>
        /// <param name="v">Second parameter</param>
        /// <returns>Value</returns>
        public static float Erlang(float y, int v)
        {
            // special cases:
            if (v == 0)
                return 1;
            if (v < 0)
                return float.NaN;

            // set:
            float t = 1, b = 1; int i;

            //series:
            for (i = 1; i < v; i++)
            {
                t *= y / i;
                b += t;
            }

            // last step and result:
            float a = t * y / i;
            return a / (a + b);
        }
        /// <summary>
        /// Returns the value of the Erlang B-function.
        /// </summary>
        /// <param name="y">First parameter</param>
        /// <param name="v">Second parameter</param>
        /// <returns>Value</returns>
        public static Complex32 Erlang(Complex32 y, int v)
        {
            // Special cases match your real version
            if (v == 0) return Complex32.One;
            if (v < 0) return Complex32.NaN;

            Complex32 t = Complex32.One;  // term for k=0
            Complex32 b = Complex32.One;  // partial sum Σ_{k=0}^{v-1}
            int i;

            // Build sum up to k = v-1:  t = y^k / k!,  b += t
            for (i = 1; i < v; i++)
            {
                t *= y / i;               // next term
                b += t;
            }

            // Last term k = v:  a = y^v / v!
            Complex32 a = t * y / i;      // here i == v
            return a / (a + b);
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
            // Only real branches supported here
            if (k != 0 && k != -1)
                return float.NaN;

            const float xm = -0.36787944117144233f; // -1/e

            // Domain checks on R
            if (k == 0)
            {
                if (x < xm) return float.NaN;   // W0 real only for x ≥ -1/e
                if (x == 0f) return 0f;
                if (x == xm) return -1f;

                // ---- ASYMPTOTIC for large x on principal branch ----
                if (x > 1000f)
                {
                    // Corless et al. expansion
                    double L1 = Math.Log(x);
                    double L2 = Math.Log(L1);
                    double invL1 = 1.0 / L1;
                    double invL1_2 = invL1 * invL1;
                    double invL1_3 = invL1_2 * invL1;

                    double v = L1 - L2
                             + L2 * invL1
                             + (L2 * (-2.0 + L2)) * 0.5 * invL1_2
                             + (L2 * (6.0 - 9.0 * L2 + 2.0 * L2 * L2)) * (1.0 / 6.0) * invL1_3;

                    return (float)v;
                }
            }
            else // k == -1
            {
                if (x < xm || x >= 0f) return float.NaN; // W_{-1} real only for -1/e ≤ x < 0
                if (x == xm) return -1f;
            }

            // Halley's method
            const float eps = 1e-8f;
            const int maxIter = 240;

            // Initial guess (simple and robust for real case)
            double w = (k == 0) ? 1f : -2f;
            // Better near the branch point x ≈ -1/e:
            if (x < -0.2f)
            {
                // w ≈ -1 ± sqrt(2 (e x + 1))
                float q = Maths.E * x + 1f;
                if (q >= 0f)
                {
                    float s = Maths.Sqrt(2f * q);
                    w = (k == 0) ? (-1f + s) : (-1f - s);
                }
            }

            for (int it = 0; it < maxIter; it++)
            {
                double e = Math.Exp(w);
                double f = w * e - x;                  // f(w) = w e^w - x
                double wp1 = w + 1f;

                // Halley denominator: e*(w+1) - (w+2)*f/(2*(w+1))
                double denom = e * wp1 - (w + 2f) * f / (2f * wp1);

                // Newton fallback if denom is unusable
                if (double.IsInfinity(denom) || denom == 0.0)
                    denom = e * wp1;

                var v = w;
                w = w - f / denom;

                if (Maths.Abs(w - v) <= eps * Maths.Abs(w))
                    break;
            }

            return (float)w;
        }
        /// <summary>
        /// Returns the value of the Lambert W-function.
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="k">Branch</param>
        /// <returns>Value</returns>
        public static Complex32 LambertW(Complex32 z, int k = 0)
        {
            // Special cases
            if (z.Real == 0f && z.Imag == 0f)
            {
                if (k == 0) return Complex32.Zero;           // W_0(0)=0
                return Complex32.NaN;                        // other branches not defined at 0
            }

            // ---- ASYMPTOTIC for large |z| on branch k ----
            // W_k(z) ≈ L1 − L2 + L2/L1 + L2(−2+L2)/(2 L1^2) + L2(6−9L2+2L2^2)/(6 L1^3)
            if (Maths.Abs(z) > 1000f)
            {
                Complex I2Pi = new Complex(0f, 2f * Maths.Pi);
                Complex L1 = Complex.Log(z) + k * I2Pi;
                Complex L2 = Complex.Log(L1);

                Complex invL1 = Complex.One / L1;
                Complex invL1_2 = invL1 * invL1;
                Complex invL1_3 = invL1_2 * invL1;

                Complex v = L1 - L2
                          + L2 * invL1
                          + (L2 * (-2f + L2)) * 0.5f * invL1_2
                          + (L2 * (6f - 9f * L2 + 2f * L2 * L2)) * (1f / 6f) * invL1_3;

                return (Complex32)v;
            }

            // Initial guess:
            // w0 ≈ L - Log(L), where L = Log(z) + i*2πk (multi-valued log).
            Complex I2Pi0 = new Complex(0f, 2f * Maths.Pi);
            Complex L = Complex.Log(z) + k * I2Pi0;
            Complex w0 = (Complex.Abs(L) < 1e-3f) ? Complex.Log(z) : (L - Complex.Log(L));
            Complex w = w0;

            float tol = 1e-8f;
            int maxIter = 240;

            // Halley's iteration for f(w)=w e^w - z
            for (int i = 0; i < maxIter; i++)
            {
                Complex ew = Complex.Exp(w);
                Complex f = w * ew - (Complex)z;
                Complex wp1 = w + Complex.One;

                // Halley denominator
                Complex denom = ew * wp1 - (w + 2f) * f / (2f * wp1);

                // Newton fallback if denom ~ 0
                if (Maths.Abs(denom) == 0f)
                    denom = ew * wp1;

                Complex wNext = w - f / denom;

                if (Maths.Abs(wNext - w) <= tol * (1f + Maths.Abs(wNext)))
                    return (Complex32)wNext;

                w = wNext;
            }

            return (Complex32)w; // last iterate
        }
        /// <summary>
        /// Returns the value of the square super-root.
        /// </summary>
        /// <param name="x">Value [1,+inf)</param>
        /// <param name="k">Branch</param>
        /// <returns>Value</returns>
        public static float Ssqrt(float x, int k = 0)
        {
            // The 2nd-order super-root, square super-root, or super square root has notation ssqrt(x).
            // It can be represented with the Lambert W-function: ssqrt(x) = log(x) / W{ log(x) }.
            float log = Maths.Log(x);
            return log / LambertW(log, k);
        }
        /// <summary>
        /// Returns the value of the square super-root.
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="k">Branch</param>
        /// <returns>Value</returns>
        public static Complex32 Ssqrt(Complex32 z, int k = 0)
        {
            // The 2nd-order super-root, square super-root, or super square root has notation ssqrt(x).
            // It can be represented with the Lambert W-function: ssqrt(x) = log(x) / W{ log(x) }.

            // z = 1 → Log(z)=0 → W_0(0)=0 → y=exp(0)=1 on principal branch
            if (z.Real == 0f && z.Imag == 0f)
                return Complex32.Zero; // y^y = 0 has solution y=0 (principal choice)

            Complex32 Lz = Maths.Log(z);  // principal log
            if (Lz.Real == 0f && Lz.Imag == 0f && k == 0)
                return Complex32.One;

            Complex32 W = LambertW(Lz, k);
            return Maths.Exp(W);
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
            if (x == 0f) return 0f;

            if (Math.Abs(x) > 6f)
            {
                double xd = x;
                double invx = 1.0 / xd;
                double x2 = xd * xd;
                double sinx2 = Math.Sin(x2);
                double cosx2 = Math.Cos(x2);

                double approx = Math.Sqrt(Math.PI / 8.0) * Math.Sign(xd)
                              + 0.5 * sinx2 * invx
                              - 0.25 * cosx2 * invx * invx * invx
                              + (3.0 / 8.0) * sinx2 * Math.Pow(invx, 5)
                              - (15.0 / 16.0) * cosx2 * Math.Pow(invx, 7);

                return (float)approx;
            }

            // term_0 = z / ((2*0)! * (4*0+1)) = z
            float eps = 1e-16f;
            int maxIter = 120;
            double s = x;
            double term = x;
            double z4 = (double)x * x; z4 *= z4; // x^4 в double

            // term_{n+1} = term_n * [ -(4n+1) z^4 / ((2n+2)(2n+1)(4n+5)) ]
            for (int n = 0; n < maxIter; n++)
            {
                double a = -(4.0 * n + 1.0);
                double b = (2.0 * n + 2.0) * (2.0 * n + 1.0) * (4.0 * n + 5.0);
                term *= a / b * z4;

                if (Math.Abs(term) < eps) break;
                s += term;
            }

            return (float)s;
        }
        /// <summary>
        /// Returns the value of the Fresnel integral C(x).
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Fresnelc(Complex32 z)
        {
            if (z.Real == 0f && z.Imag == 0f) return Complex32.Zero;
            if (z.Imag == 0f) return new Complex32(Fresnelc(z.Real), 0f);

            double az = Math.Sqrt((double)z.Real * z.Real + (double)z.Imag * z.Imag);

            if (az > 6.0)
            {
                double sgnInf = Math.Sign((double)z.Real + (double)z.Imag);

                var Z = new Complex(z.Real, z.Imag);
                var invZ = Complex.One / Z;
                var Z2 = Z * Z;
                var sinZ2 = Complex.Sin(Z2);
                var cosZ2 = Complex.Cos(Z2);

                var approx = Math.Sqrt(Math.PI / 8.0) * sgnInf
                           + 0.5 * sinZ2 * invZ
                           - 0.25 * cosZ2 * invZ * invZ * invZ
                           + (3.0 / 8.0) * sinZ2 * invZ * invZ * invZ * invZ * invZ
                           - (15.0 / 16.0) * cosZ2 * invZ * invZ * invZ * invZ * invZ * invZ * invZ;

                return new Complex32((float)approx.Real, (float)approx.Imaginary);
            }

            // term_0 = z / ((2*0)! * (4*0+1)) = z
            float eps = 1e-16f;
            int maxIter = 120;
            var s = new Complex(z.Real, z.Imag);
            var term = s;
            var z4 = s * s * s * s;

            // term_{n+1} = term_n * [ -(4n+1) z^4 / ((2n+2)(2n+1)(4n+5)) ]
            for (int n = 0; n < maxIter; n++)
            {
                double a = -(4.0 * n + 1.0);
                double b = (2.0 * n + 2.0) * (2.0 * n + 1.0) * (4.0 * n + 5.0);
                term *= a / b * z4;

                if (Complex.Abs(term) < eps) break;
                s += term;
            }

            return new Complex32((float)s.Real, (float)s.Imaginary);
        }
        /// <summary>
        /// Returns the value of the Fresnel integral S(x).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Fresnels(float x)
        {
            if (x == 0f) return 0f;

            if (Math.Abs(x) > 6f)
            {
                double xd = x;
                double invx = 1.0 / xd;
                double x2 = xd * xd;
                double sinx2 = Math.Sin(x2);
                double cosx2 = Math.Cos(x2);

                double approx = Math.Sqrt(Math.PI / 8.0) * Math.Sign(xd)
                              - 0.5 * cosx2 * invx
                              - 0.25 * sinx2 * invx * invx * invx
                              - (3.0 / 8.0) * cosx2 * Math.Pow(invx, 5)
                              - (15.0 / 16.0) * sinx2 * Math.Pow(invx, 7);

                return (float)approx;
            }

            // term_0 = z^3 / 3 = z^{4*0+3} / ((2*0+1)! (4*0+3))
            float eps = 1e-16f;
            int maxIter = 120;
            double z2 = (double)x * x;
            double s = z2 * x / 3.0;
            double term = s;
            double z4 = z2 * z2;

            // term_{n+1} = term_n * [ -(4n+3) z^4 / ((2n+3)(2n+2)(4n+7)) ]
            for (int n = 0; n < maxIter; n++)
            {
                double a = -(4.0 * n + 3.0);
                double b = (2.0 * n + 3.0) * (2.0 * n + 2.0) * (4.0 * n + 7.0);
                term *= (a / b) * z4;

                if (Math.Abs(term) < eps) break;
                s += term;
            }

            return (float)s;
        }
        /// <summary>
        /// Returns the value of the Fresnel integral S(x).
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Fresnels(Complex32 z)
        {
            if (z.Real == 0f && z.Imag == 0f) return Complex32.Zero;
            if (z.Imag == 0f) return new Complex32(Fresnels(z.Real), 0f);

            double az = Math.Sqrt((double)z.Real * z.Real + (double)z.Imag * z.Imag);

            if (az > 6.0)
            {
                double sgnInf = Math.Sign((double)z.Real + (double)z.Imag);

                var Z = new Complex(z.Real, z.Imag);
                var invZ = Complex.One / Z;
                var Z2 = Z * Z;
                var sinZ2 = Complex.Sin(Z2);
                var cosZ2 = Complex.Cos(Z2);

                var approx = Math.Sqrt(Math.PI / 8.0) * sgnInf
                           - 0.5 * cosZ2 * invZ
                           - 0.25 * sinZ2 * invZ * invZ * invZ
                           - (3.0 / 8.0) * cosZ2 * invZ * invZ * invZ * invZ * invZ
                           - (15.0 / 16.0) * sinZ2 * invZ * invZ * invZ * invZ * invZ * invZ * invZ;

                return new Complex32((float)approx.Real, (float)approx.Imaginary);
            }

            // term_0 = z^3 / 3 = z^{4*0+3} / ((2*0+1)! (4*0+3))
            float eps = 1e-16f;
            int maxIter = 120;
            var Zc = new Complex(z.Real, z.Imag);
            var z2 = Zc * Zc;
            var s = z2 * Zc / 3.0;
            var term = s;
            var z4 = z2 * z2;

            // term_{n+1} = term_n * [ -(4n+3) z^4 / ((2n+3)(2n+2)(4n+7)) ]
            for (int n = 0; n < maxIter; n++)
            {
                double a = -(4.0 * n + 3.0);
                double b = (2.0 * n + 3.0) * (2.0 * n + 2.0) * (4.0 * n + 7.0);
                term *= a / b * z4;

                if (Complex.Abs(term) < eps) break;
                s += term;
            }

            return new Complex32((float)s.Real, (float)s.Imaginary);
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
            if (float.IsNaN(h) || float.IsNaN(a)) return float.NaN;
            if (a == 0f) return 0f;
            if (h == 0f) return INV_2PI * Maths.Atan(a);

            float sign = a >= 0f ? 1f : -1f;
            float L = Math.Abs(a);

            int n = 1024;
            float hstep = L / n;

            float f(float t)
            {
                float t2 = t * t;
                return Maths.Exp(-0.5f * h * h * (1f + t2)) / (1f + t2);
            }

            float sum = f(0f) + f(L);
            float s4 = 0f, s2 = 0f;
            for (int i = 1; i < n; i++)
            {
                float ti = i * hstep;
                if ((i & 1) == 1) s4 += f(ti); else s2 += f(ti);
            }
            float integral = hstep / 3f * (sum + 4f * s4 + 2f * s2);
            return sign * INV_2PI * integral;
        }
        /// <summary>
        /// Returns the value of the Owen T function.
        /// </summary>
        /// <param name="h">First value</param>
        /// <param name="a">Second value</param>
        /// <returns>Value</returns>
        public static Complex32 Owen(Complex32 h, Complex32 a)
        {
            if (float.IsNaN(h.Real) || float.IsNaN(h.Imag) ||
                float.IsNaN(a.Real) || float.IsNaN(a.Imag))
                return Complex32.NaN;

            if (a.Real == 0f && a.Imag == 0f) return Complex32.Zero;
            if (h.Real == 0f && h.Imag == 0f)
            {
                return INV_2PI * Maths.Atan(a);
            }

            int n = 1024;
            float ds = 1f / n;

            Complex32 F(float s)
            {
                Complex32 sa = s * a;
                Complex32 denom = Complex32.One + sa * sa;
                Complex32 expo = Maths.Exp(-0.5f * h * h * denom);
                return a * (expo / denom); // dt = a ds
            }

            Complex32 sum = F(0f) + F(1f);
            Complex32 s4 = Complex32.Zero;
            Complex32 s2 = Complex32.Zero;
            for (int i = 1; i < n; i++)
            {
                float si = i * ds;
                if ((i & 1) == 1) s4 += F(si); else s2 += F(si);
            }
            Complex32 integral = ds / 3f * (sum + 4f * s4 + 2f * s2);
            return INV_2PI * integral;
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
            if (float.IsNaN(s)) return float.NaN;
            if (s == 1f) return float.PositiveInfinity;     // simple pole

            // Re(s) > 0: use eta(s)/(1-2^{1-s})
            if (s > 0f)
                return ZetaHasse(s);

            // Functional equation for s ≤ 0
            float twoPowS = Maths.Exp(s * Maths.Log(2f));
            float piPow = Maths.Exp((s - 1f) * Maths.Log(Maths.Pi));
            float sinTerm = Maths.Sin(0.5f * Maths.Pi * s);
            float gamma = Special.Gamma(1f - s);

            float zeta1ms = ZetaHasse(1f - s);          // now 1-s ≥ 1
            return twoPowS * piPow * sinTerm * gamma * zeta1ms;
        }
        /// <summary>
        /// Returns the value of the Riemann zeta ζ(s) on the principal branch (complex s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Zeta(Complex32 s)
        {
            // explicit pole at s = 1 + 0i
            if (s.Imag == 0f && s.Real == 1f)
                return Complex32.NaN;

            if (s.Real > 0f && !(s.Imag == 0f && s.Real == 1f))
                return ZetaHasse(s);

            // Functional equation
            Complex32 two = new Complex32(2f, 0f);
            Complex32 pi = new Complex32(Maths.Pi, 0f);

            Complex32 twoPowS = Maths.Exp(s * Maths.Log(two));
            Complex32 piPow = Maths.Exp((s - Complex32.One) * Maths.Log(pi));
            Complex32 sinTerm = Maths.Sin(0.5f * pi * s);
            Complex32 gamma = Special.Gamma(Complex32.One - s);

            Complex32 zeta1ms = ZetaHasse(Complex32.One - s);
            return twoPowS * piPow * sinTerm * gamma * zeta1ms;
        }

        #region Private methods (helpers)

        // ---------- helpers (Dirichlet eta series, principal branch) ----------

        /// <summary>
        /// Evaluates the Dirichlet eta function for real inputs via an alternating series.
        /// </summary>
        /// <param name="s">Function argument</param>
        [Obsolete]
        private static float ZetaEta(float s)
        {
            // ζ(s) = η(s) / (1 - 2^{1-s}), valid for s>0, s≠1
            float denom = 1f - Maths.Exp((1f - s) * Maths.Log(2f));
            // Handle s≈1 numerically? Here API has exact s; s==1 already excluded.

            const float eps = 1e-16f;
            const int maxIter = 100000;

            double sum = 0.0; // accumulate in double for a bit more headroom
            for (int n = 1; n <= maxIter; n++)
            {
                double term = 1.0 / Math.Pow(n, s);
                if ((n & 1) == 0) term = -term; // alternating
                sum += term;
                if (Math.Abs(term) < eps) break;
            }
            return (float)(sum / denom);
        }

        /// <summary>
        /// Evaluates the Dirichlet eta function for complex inputs using an alternating series.
        /// </summary>
        /// <param name="s">Function argument</param>
        [Obsolete]
        private static Complex32 ZetaEta(Complex32 s)
        {
            // ζ(s) = η(s) / (1 - 2^{1-s}), valid for Re(s)>0, s≠1
            Complex32 two = new Complex32(2f, 0f);
            Complex32 denom = Complex32.One - Maths.Exp((Complex32.One - s) * Maths.Log(two));

            const float eps = 1e-16f;
            const int maxIter = 10000000;

            Complex32 sum = Complex32.Zero;
            for (int n = 1; n <= maxIter; n++)
            {
                // term = (-1)^{n-1} / n^s  = (-1)^{n-1} * exp( -s * log n )
                Complex32 term = Maths.Exp(-s * new Complex32(Maths.Log(n), 0f));
                if ((n & 1) == 0) term = -term;
                sum += term;
                if (Maths.Abs(term) < eps) break;
            }
            return sum / denom;
        }

        /// <summary>
        /// Computes the Riemann zeta function for real inputs using Hasse's series.
        /// </summary>
        /// <param name="s">Function argument</param>
        private static float ZetaHasse(float s)
        {
            // ζ(s) = S(s) / (1 - 2^{1-s}), with S(s) the Hasse inner sum
            float denom = 1f - Maths.Exp((1f - s) * Maths.Log(2f));
            // handle near-pole numerically? caller excludes s=1 exactly.

            double S = 0.0;
            float eps = 1e-16f;
            int maxK = 512;

            for (int k = 0; k < maxK; k++)
            {
                // inner sum: A_k(s) = Σ_{q=0..k} (-1)^q C(k,q) (q+1)^{-s}
                double Ak = 0.0;
                double C = 1.0; // C(k,0)

                for (int q = 0; q <= k; q++)
                {
                    double term = C * Math.Pow(q + 1.0, -s);
                    if ((q & 1) == 1) term = -term;
                    Ak += term;

                    // next binomial C(k,q+1)
                    C *= (k - q) / (q + 1.0);
                }

                double Tk = Ak / Math.Pow(2.0, k + 1.0);
                S += Tk;

                if (Math.Abs(Tk) < eps * (1.0 + Math.Abs(S))) break;
            }
            return (float)(S / denom);
        }

        /// <summary>
        /// Computes the Riemann zeta function for complex inputs using Hasse's series.
        /// </summary>
        /// <param name="s">Function argument</param>
        private static Complex32 ZetaHasse(Complex32 s)
        {
            Complex32 two = new Complex32(2f, 0f);
            Complex32 denom = Complex32.One - Maths.Exp((Complex32.One - s) * Maths.Log(two));

            Complex32 S = Complex32.Zero;
            float eps = 1e-16f;
            int maxK = 512;

            for (int k = 0; k < maxK; k++)
            {
                Complex32 Ak = Complex32.Zero;
                // binomial C(k,q) iteratively
                Complex32 C = Complex32.One; // real-positive actually, but keep as complex-friendly
                for (int q = 0; q <= k; q++)
                {
                    // (q+1)^{-s} = exp( -s * log(q+1) )
                    Complex32 pow = Maths.Exp(-s * new Complex32(Maths.Log(q + 1), 0f));
                    Complex32 term = C * pow;
                    if ((q & 1) == 1) term = -term;
                    Ak += term;

                    // C(k,q+1) = C(k,q) * (k - q)/(q + 1)
                    float ratio = (q < k) ? (float)(k - q) / (q + 1f) : 0f;
                    C *= new Complex32(ratio, 0f);
                }

                Complex32 Tk = Ak / Maths.Exp(new Complex32((k + 1f) * Maths.Log(2f), 0f));
                S += Tk;

                if (Maths.Abs(Tk) < eps * (1f + Maths.Abs(S))) break;
            }
            return S / denom;
        }

        #endregion

        #endregion

        #region Gamma functions

        /// <summary>
        /// Returns the value of the Euler Gamma function: Г(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Gamma(float x)
        {
            if (float.IsNaN(x)) return float.NaN;

            // Poles on the real axis
            if (x <= 0f && x == Maths.Round(x)) // integer check
            {
                int xi = (int)Maths.Round(x);
                if (xi <= 0) return float.NaN;
            }

            // Reflection for better accuracy and negative non-integers
            if (x < 0.5f)
            {
                float sinpix = Maths.Sin(Maths.Pi * x);
                if (sinpix == 0f) return float.NaN; // pole
                return Maths.Pi / (sinpix * Gamma(1f - x));
            }

            return (float)GammaLanczos((double)x);
        }
        /// <summary>
        /// Returns the value of the Euler Gamma function: Г(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Gamma(Complex32 x)
        {
            // Detect real negative integer poles exactly
            if (x.Imag == 0f)
            {
                float xr = x.Real;
                if (xr <= 0f && xr == Maths.Round(xr)) return Complex32.NaN;
            }

            if (x.Real < 0.5f)
            {
                Complex32 sinpiz = Maths.Sin(new Complex32(Maths.Pi, 0f) * x);
                if (Maths.Abs(sinpiz) == 0f) return Complex32.NaN;
                return new Complex32(Maths.Pi, 0f) / (sinpiz * Gamma(Complex32.One - x));
            }

            return GammaLanczos(x);
        }

        /// <summary>
        /// Returns the value of the natural logarithm of the Euler Gamma function: ln[Г(z)].
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float LogGamma(float x)
        {
            if (float.IsNaN(x)) return float.NaN;

            // Poles on the real axis
            if (x <= 0f && x == Maths.Round(x)) return float.NaN;

            if (x < 0.5f)
            {
                // log|Γ(x)| = log(π) - log|sin(πx)| - log|Γ(1-x)|
                float sinpix = Maths.Sin(Maths.Pi * x);
                if (sinpix == 0f) return float.NaN;
                return Maths.Log(Maths.Pi) - Maths.Log(Maths.Abs(sinpix)) - LogGamma(1f - x);
            }

            return (float)LogGammaLanczos((double)x);
        }
        /// <summary>
        /// Returns the value of the natural logarithm of the Euler Gamma function: ln[Г(z)].
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 LogGamma(Complex32 z)
        {
            // Real negative integers: simple poles
            if (z.Imag == 0f)
            {
                float xr = z.Real;
                if (xr <= 0f && xr == Maths.Round(xr))
                    return Complex32.NaN;
            }

            if (z.Real < 0.5f)
            {
                // log Γ(z) = log π − log sin(πz) − log Γ(1−z)  (principal branches)
                Complex32 sinpiz = Maths.Sin(new Complex32(Maths.Pi, 0f) * z);
                if (Maths.Abs(sinpiz) == 0f) return Complex32.NaN;
                return Maths.Log(new Complex32(Maths.Pi, 0f)) - Maths.Log(sinpiz) - LogGamma(Complex32.One - z);
            }

            return LogGammaLanczos(z);
        }

        /// <summary>
        /// Returns the value of the Digamma function: ψ(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float DiGamma(float x)
        {
            if (float.IsNaN(x)) return float.NaN;

            // Real poles at non-positive integers
            if (x <= 0f && x == Maths.Round(x)) return float.NaN;

            // Reflection to avoid small/negative arguments
            if (x < 0.5f)
            {
                float pix = Maths.Pi * x;
                float sin = Maths.Sin(pix);
                if (sin == 0f) return float.NaN;             // pole
                                                             // ψ(x) = ψ(1-x) - π cot(πx)
                return DiGamma(1f - x) - Maths.Pi * (Maths.Cos(pix) / sin);
            }

            // Shift up to a safe region for the asymptotic
            float acc = 0f;
            float z = x;
            while (z < 8f)
            {
                acc -= 1f / z;
                z += 1f;
            }

            // Bernoulli asymptotic: ψ(z) ~ ln z - 1/(2z) - 1/(12 z^2) + 1/(120 z^4) - 1/(252 z^6) + 1/(240 z^8) - 1/(132 z^10)
            float inv = 1f / z;
            float inv2 = inv * inv;
            float res = Maths.Log(z) - 0.5f * inv
                        - (1f / 12f) * inv2
                        + (1f / 120f) * (inv2 * inv2)
                        - (1f / 252f) * (inv2 * inv2 * inv2)
                        + (1f / 240f) * (inv2 * inv2 * inv2 * inv2)
                        - (1f / 132f) * (inv2 * inv2 * inv2 * inv2 * inv2);

            return res + acc;
        }
        /// <summary>
        /// Returns the value of the Digamma function: ψ(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 DiGamma(Complex32 x)
        {
            // Real negative integers: poles
            if (x.Imag == 0f)
            {
                float xr = x.Real;
                if (xr <= 0f && xr == Maths.Round(xr))
                    return Complex32.NaN;
            }

            Complex32 z = x;

            // Reflection for better behavior near the poles/left half-plane
            if (z.Real < 0.5f)
            {
                Complex32 piz = new Complex32(Maths.Pi, 0f) * z;
                Complex32 sin = Maths.Sin(piz);
                if (Maths.Abs(sin) == 0f) return Complex32.NaN;
                // ψ(z) = ψ(1 - z) - π cot(π z) = ψ(1 - z) - π * cos(π z)/sin(π z)
                Complex32 cot = Maths.Cos(piz) / sin;
                return DiGamma(Complex32.One - z) - new Complex32(Maths.Pi, 0f) * cot;
            }

            // Shift up to Re(z) ≥ 8
            Complex32 acc = Complex32.Zero;
            while (z.Real < 8f)
            {
                acc -= Complex32.One / z;
                z += Complex32.One;
            }

            // Bernoulli asymptotic in the complex plane
            Complex32 inv = Complex32.One / z;
            Complex32 inv2 = inv * inv;

            Complex32 res = Maths.Log(z) - 0.5f * inv
                            - (1f / 12f) * inv2
                            + (1f / 120f) * (inv2 * inv2)
                            - (1f / 252f) * (inv2 * inv2 * inv2)
                            + (1f / 240f) * (inv2 * inv2 * inv2 * inv2)
                            - (1f / 132f) * (inv2 * inv2 * inv2 * inv2 * inv2);

            return res + acc;
        }

        /// <summary>
        /// Returns the value of the Trigamma function: ψ1(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float TriGamma(float x)
        {
            if (float.IsNaN(x)) return float.NaN;

            // Real poles at non-positive integers
            if (x <= 0f && x == Maths.Round(x)) return float.NaN;

            // Reflection for better behavior in left half-plane
            if (x < 0.5f)
            {
                float pix = Maths.Pi * x;
                float sin = Maths.Sin(pix);
                if (sin == 0f) return float.NaN;
                float csc2 = 1f / (sin * sin);
                return (Maths.Pi * Maths.Pi) * csc2 - TriGamma(1f - x);
            }

            // Shift up to a safe region for asymptotics
            float z = x;
            float acc = 0f;
            while (z < 8f)
            {
                acc += 1f / (z * z);
                z += 1f;
            }

            // Bernoulli asymptotic
            float inv = 1f / z;
            float inv2 = inv * inv;
            float inv3 = inv2 * inv;
            float inv5 = inv3 * inv2;
            float inv7 = inv5 * inv2;
            float inv9 = inv7 * inv2;
            float inv11 = inv9 * inv2;

            float res = inv
                       + 0.5f * inv2
                       + (1f / 6f) * inv3
                       - (1f / 30f) * inv5
                       + (1f / 42f) * inv7
                       - (1f / 30f) * inv9
                       + (5f / 66f) * inv11;

            return res + acc;
        }
        /// <summary>
        /// Returns the value of the Trigamma function: ψ1(z).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 TriGamma(Complex32 x)
        {
            // Detect real negative integer poles exactly
            if (x.Imag == 0f)
            {
                float xr = x.Real;
                if (xr <= 0f && xr == Maths.Round(xr))
                    return Complex32.NaN;
            }

            Complex32 z = x;

            // Reflection
            if (z.Real < 0.5f)
            {
                Complex32 piz = new Complex32(Maths.Pi, 0f) * z;
                Complex32 sin = Maths.Sin(piz);
                if (Maths.Abs(sin) == 0f) return Complex32.NaN;

                Complex32 csc2 = Complex32.One / (sin * sin);
                return new Complex32(Maths.Pi * Maths.Pi, 0f) * csc2 - TriGamma(Complex32.One - z);
            }

            // Shift up to Re(z) ≥ 8
            Complex32 acc = Complex32.Zero;
            while (z.Real < 8f)
            {
                acc += Complex32.One / (z * z);
                z += Complex32.One;
            }

            // Bernoulli asymptotic
            Complex32 inv = Complex32.One / z;
            Complex32 inv2 = inv * inv;
            Complex32 inv3 = inv2 * inv;
            Complex32 inv5 = inv3 * inv2;
            Complex32 inv7 = inv5 * inv2;
            Complex32 inv9 = inv7 * inv2;
            Complex32 inv11 = inv9 * inv2;

            Complex32 res = inv
                           + 0.5f * inv2
                           + (1f / 6f) * inv3
                           - (1f / 30f) * inv5
                           + (1f / 42f) * inv7
                           - (1f / 30f) * inv9
                           + (5f / 66f) * inv11;

            return res + acc;
        }

        /// <summary>
        /// Returns the value of the incomplete upper Gamma function: Q(s, x) = Γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float GammaQ(float s, float x)
        {
            if (float.IsNaN(s) || float.IsNaN(x)) return float.NaN;
            if (s <= 0f && s == Maths.Round(s)) return float.NaN;    // poles of Γ(s)
            if (x < 0f) return float.NaN;

            if (x == 0f) return 1f;          // Q(s,0)=1 (for s>0)
            if (float.IsPositiveInfinity(x)) return 0f;

            // choose method
            if (x < s + 1f)
            {
                // series for P(s,x), then Q = 1 - P
                float P = LowerRegGammaSeries(s, x);
                return 1f - P;
            }
            else
            {
                // continued fraction for Q(s,x)
                float logPref = s * Maths.Log(x) - x - Special.LogGamma(s); // log( e^{-x} x^s / Γ(s) )
                float h = UpperGammaCF(s, x);                 // CF value
                return Maths.Exp(logPref) * h;
            }
        }
        /// <summary>
        /// Returns the value of the incomplete upper Gamma function: Q(s, x) = Γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 GammaQ(Complex32 s, Complex32 x)
        {
            // Poles at real non-positive integers s
            if (x.Real == 0f && x.Imag == 0f) return Complex32.One; // Q(s,0)=1 (analytic continuation)
            if (float.IsNaN(s.Real) || float.IsNaN(s.Imag) ||
                float.IsNaN(x.Real) || float.IsNaN(x.Imag))
                return Complex32.NaN;

            // heuristic split like real case
            if (Maths.Abs(x) < Maths.Abs(s) + 1f)
            {
                Complex32 P = LowerRegGammaSeries(s, x);
                return Complex32.One - P;
            }
            else
            {
                Complex32 logPref = s * Maths.Log(x) - x - Special.LogGamma(s);
                Complex32 h = UpperGammaCF(s, x);
                return Maths.Exp(logPref) * h;
            }
        }

        /// <summary>
        /// Returns the value of an incomplete lower Gamma function: P(s, x) = γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float GammaP(float s, float x)
        {
            if (float.IsNaN(s) || float.IsNaN(x)) return float.NaN;
            if (x < 0f) return float.NaN;
            // Poles of Γ(s) on the real line:
            if (s <= 0f && s == Maths.Round(s)) return float.NaN;

            if (x == 0f) return 0f;                        // P(s,0) = 0  (s>0)
            if (float.IsPositiveInfinity(x)) return 1f;    // P(s,∞) = 1

            // Choose method:
            if (x < s + 1f)
            {
                // Direct series for P(s,x)
                return LowerRegGammaSeries(s, x);
            }
            else
            {
                // Continued fraction for Q(s,x), then P = 1 - Q
                float logPref = s * Maths.Log(x) - x - Special.LogGamma(s); // log(e^{-x} x^s / Γ(s))
                float cf = UpperGammaCF(s, x);                // CF approximates Γ(s,x)/(e^{-x} x^s)
                float Q = Maths.Exp(logPref) * cf;
                return 1f - Q;
            }
        }
        /// <summary>
        /// Returns the value of an incomplete lower Gamma function: P(s, x) = γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 GammaP(Complex32 s, Complex32 x)
        {
            // Poles of Γ(s) at real non-positive integers:
            if (x.Real == 0f && x.Imag == 0f) return Complex32.Zero; // P(s,0)=0 (analytic continuation)
            if ((s.Imag == 0f) && s.Real <= 0f && s.Real == Maths.Round(s.Real))
                return Complex32.NaN;

            // Heuristic split mirroring the real case
            if (Maths.Abs(x) < Maths.Abs(s) + 1f)
            {
                return LowerRegGammaSeries(s, x);
            }
            else
            {
                Complex32 logPref = s * Maths.Log(x) - x - Special.LogGamma(s);
                Complex32 cf = UpperGammaCF(s, x);
                Complex32 Q = Maths.Exp(logPref) * cf;
                return Complex32.One - Q;
            }
        }

        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float GammaIncomplete(float s, float x)
        {
            if (float.IsNaN(s) || float.IsNaN(x)) return float.NaN;
            if (x < 0f) return float.NaN;

            // Poles of Γ(s) on the real axis (we restrict to s > 0 in this real overload)
            if (s <= 0f && s == Maths.Round(s)) return float.NaN;

            if (x == 0f) return 0f;                       // γ(s,0) = 0  (s>0)
            if (float.IsPositiveInfinity(x)) return Special.Gamma(s); // γ(s,∞) = Γ(s)

            if (x < s + 1f)
            {
                // Series for γ(s,x) via P(s,x) series, then multiply by Γ(s)
                float P = LowerRegGammaSeries(s, x);    // P(s,x)
                return Special.Gamma(s) * P;                          // γ = Γ * P
            }
            else
            {
                // γ(s,x) = Γ(s) − Γ(s,x),  Γ(s,x) ≈ e^{−x} x^s * CF
                float logPref = s * Maths.Log(x) - x;                  // log(e^{−x} x^s)
                float cf = UpperGammaCF(s, x);           // CF ≈ Γ(s,x)/(e^{−x} x^s)
                float upper = Maths.Exp(logPref) * cf;
                return Special.Gamma(s) - upper;
            }
        }
        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 GammaIncomplete(Complex32 s, Complex32 x)
        {
            if ((float.IsNaN(s.Real) || float.IsNaN(s.Imag)) ||
                (float.IsNaN(x.Real) || float.IsNaN(x.Imag)))
                return Complex32.NaN;

            if (x.Real == 0f && x.Imag == 0f) return Complex32.Zero; // γ(s,0) = 0

            // Poles of Γ(s) on the principal branch (real negative integers)
            if (s.Imag == 0f)
            {
                float sr = s.Real;
                if (sr <= 0f && sr == Maths.Round(sr))
                    return Complex32.NaN;
            }

            // Heuristic split mirroring the real case
            if (Maths.Abs(x) < Maths.Abs(s) + 1f)
            {
                // γ = Γ * P, with P via series
                Complex32 P = LowerRegGammaSeries(s, x);
                return Special.Gamma(s) * P;
            }
            else
            {
                // γ = Γ − Γ(s,x), with Γ(s,x) ≈ e^{−x} x^s * CF
                Complex32 logPref = s * Maths.Log(x) - x;         // log(e^{−x} x^s)
                Complex32 cf = UpperGammaCF(s, x);
                Complex32 upper = Maths.Exp(logPref) * cf;
                return Special.Gamma(s) - upper;
            }
        }

        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x) (complemented).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float GammaIncompleteComplemented(float s, float x)
        {
            // Γ(s,x) = Γ(s) - γ(s,x)
            return Special.Gamma(s) - Special.GammaIncomplete(s, x);
        }
        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x) (complemented).
        /// </summary>
        /// <param name="s">Value</param>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 GammaIncompleteComplemented(Complex32 s, Complex32 x)
        {
            return Special.Gamma(s) - Special.GammaIncomplete(s, x);
        }

        #region Private methods (helpers)

        // ---------------------- helpers ----------------------

        private const float G = 7f;
        private const float SQRT_TWO_PI = 2.5066282746310005024f; // √(2π)

        // Wikipedia/Numerical Recipes coefficients (cast to float)
        private static readonly float[] LANCZOS = new float[]
        {
            0.99999999999980993f,
            676.5203681218851f,
            -1259.1392167224028f,
            771.32342877765313f,
            -176.61502916214059f,
            12.507343278686905f,
            -0.13857109526572012f,
            0.0000099843695780195716f,
            0.00000015056327351493116f
        };

        /// <summary>
        /// Lanczos approximation (real) valid for Re(z) > 0.5. Computes principal branch.
        /// </summary>
        private static double GammaLanczos(double x)
        {
            // Γ(x) ≈ √(2π) * (t)^{x-0.5} * e^{-t} * A(x), where t = x + g - 0.5, A(x)=a0 + Σ a_k/(x-1+k)
            double z = x - 1.0;
            double sum = LANCZOS[0];
            for (int k = 1; k < LANCZOS.Length; k++)
                sum += LANCZOS[k] / (z + k);

            double t = z + G + 0.5;
            // Use log form to reduce overflow/underflow
            double logGamma = Math.Log(SQRT_TWO_PI * sum) + (z + 0.5) * Math.Log(t) - t;
            return Math.Exp(logGamma);
        }

        /// <summary>
        /// Lanczos approximation (complex) valid for Re(z) > 0.5. Computes principal branch.
        /// </summary>
        private static Complex32 GammaLanczos(Complex32 x)
        {
            Complex32 z = x - Complex32.One;
            Complex32 sum = new Complex32(LANCZOS[0], 0f);
            for (int k = 1; k < LANCZOS.Length; k++)
                sum += LANCZOS[k] / (z + new Complex32(k, 0f));

            Complex32 t = z + new Complex32(G + 0.5f, 0f);
            Complex32 logPart = Maths.Log(new Complex32(SQRT_TWO_PI, 0f) * sum)
                              + (z + new Complex32(0.5f, 0f)) * Maths.Log(t)
                              - t;
            return Maths.Exp(logPart);
        }

        /// <summary>
        /// Lanczos log-gamma for real x with Re(x) ≥ 0.5 (principal branch).
        /// </summary>
        private static double LogGammaLanczos(double x)
        {
            double z = x - 1.0;
            double sum = LANCZOS[0];
            for (int k = 1; k < LANCZOS.Length; k++)
                sum += LANCZOS[k] / (z + k);

            double t = z + G + 0.5;
            // log Γ(x) ≈ log(√(2π)*sum) + (z+0.5) log t − t
            return Math.Log(SQRT_TWO_PI * sum) + (z + 0.5) * Math.Log(t) - t;
        }

        /// <summary>
        /// Lanczos log-gamma for complex z with Re(z) ≥ 0.5 (principal branch).
        /// </summary>
        private static Complex32 LogGammaLanczos(Complex32 x)
        {
            Complex32 z = x - Complex32.One;
            Complex32 sum = new Complex32(LANCZOS[0], 0f);
            for (int k = 1; k < LANCZOS.Length; k++)
                sum += LANCZOS[k] / (z + new Complex32(k, 0f));

            Complex32 t = z + new Complex32(G + 0.5f, 0f);
            // log Γ(x) ≈ log(√(2π)*sum) + (z+0.5) log t − t
            return Maths.Log(new Complex32(SQRT_TWO_PI, 0f) * sum)
                 + (z + new Complex32(0.5f, 0f)) * Maths.Log(t)
                 - t;
        }

        // ---------------- helpers: series for P(s,x) ----------------
        // P(s,x) = e^{-x} x^{s} * Σ_{n>=0} [ x^n / (s(s+1)...(s+n)) ],
        // with t0 = 1/s and t_{n+1} = t_n * x/(s+n+1).
        /// <summary>
        /// Series expansion for the lower regularized incomplete gamma function.
        /// </summary>
        /// <param name="s">Shape parameter</param>
        /// <param name="x">Upper limit of integration</param>
        private static float LowerRegGammaSeries(float s, float x)
        {
            float term = 1f / s;    // n=0
            float sum = term;
            float eps = 1e-16f;
            int maxIter = 200;

            for (int n = 0; n < maxIter; n++)
            {
                term *= x / (s + n + 1f);
                sum += term;
                if (Maths.Abs(term) < eps * (1f + Maths.Abs(sum))) break;
            }
            // e^{-x} x^{s}
            float pref = Maths.Exp(-x + s * Maths.Log(x));
            return pref * sum;
        }

        /// <summary>
        /// Series expansion for the lower regularized incomplete gamma function with complex parameters.
        /// </summary>
        /// <param name="s">Shape parameter</param>
        /// <param name="x">Upper limit of integration</param>
        private static Complex32 LowerRegGammaSeries(Complex32 s, Complex32 x)
        {
            Complex32 term = Complex32.One / s; // n=0
            Complex32 sum = term;
            float eps = 1e-16f;
            int maxIter = 200;

            for (int n = 0; n < maxIter; n++)
            {
                term *= x / (s + (n + 1f));
                sum += term;
                if (Maths.Abs(term) < eps * (1f + Maths.Abs(sum))) break;
            }
            Complex32 pref = Maths.Exp(-x + s * Maths.Log(x));
            return pref * sum;
        }

        // --------------- helpers: continued fraction for Q(s,x) ---------------
        // Lentz’s algorithm for the continued fraction representation of Γ(s,x):
        // Q(s,x) = e^{-x} x^{s} / Γ(s) * h, where
        //   b0 = x + 1 - s,  h = 1/b0,
        //   for i=1..: a_i = i*(s - i),  b_i = b_{i-1} + 2,
        //   update via Lentz with c,d and delta = c*d.
        /// <summary>
        /// Continued fraction evaluation for the upper incomplete gamma function.
        /// </summary>
        /// <param name="s">Shape parameter</param>
        /// <param name="x">Lower limit of integration</param>
        private static float UpperGammaCF(float s, float x)
        {
            const float tiny = 1e-30f;

            float b = x + 1f - s;
            float c = 1f / tiny;
            float d = 1f / b;
            float h = d;

            float eps = 1e-16f;
            int maxIter = 200;

            for (int i = 1; i <= maxIter; i++)
            {
                float a = i * (s - i);      // a_i
                                            // d = 1 / (b + a*d),  c = b + a/c
                d = 1f / (b + a * d);
                c = b + a / c;
                float delta = c * d;
                h *= delta;

                b += 2f;
                if (Maths.Abs(delta - 1f) < eps) break;
            }
            return h;
        }

        /// <summary>
        /// Continued fraction evaluation for the upper incomplete gamma function with complex parameters.
        /// </summary>
        /// <param name="s">Shape parameter</param>
        /// <param name="x">Lower limit of integration</param>
        private static Complex32 UpperGammaCF(Complex32 s, Complex32 x)
        {
            const float tiny = 1e-30f;

            Complex32 b = x + Complex32.One - s;
            Complex32 c = new Complex32(1f / tiny, 0f);
            Complex32 d = Complex32.One / b;
            Complex32 h = d;

            float eps = 1e-16f;
            int maxIter = 200;

            for (int i = 1; i <= maxIter; i++)
            {
                Complex32 ii = new Complex32(i, 0f);
                Complex32 a = ii * (s - ii); // a_i

                d = Complex32.One / (b + a * d);
                c = b + a / c;
                Complex32 delta = c * d;
                h *= delta;

                b += 2f; // b increases by 2 each step
                if (Maths.Abs(delta - Complex32.One) < eps) break;
            }
            return h;
        }

        #endregion

        #endregion

        #region Generalized error function
        /// <summary>
        /// Returns the value of the generalized error function.
        /// </summary>
        /// <param name="x">Value (0, +inf)</param>
        /// <param name="n">Order [0, +inf)</param>
        /// <returns>Value</returns>
        public static float Gerf(float x, int n)
        {
            // Generalized error functions:
            if (n < 0)
                return float.NaN; // singular values

            else if (n == 0)
                return x / Maths.E / SQRT_PI; // E0(x)

            else if (n == 1)
                return (1 - Maths.Exp(-x)) / SQRT_PI; // E1(x)

            else if (n == 2)
                return Erf(x); // E2(x) = erf(x)

            // En(x) for all x > 0
            if (x > 0)
            {
                float p = 1.0f / n;
                float w = 1.0f / SQRT_PI;
                float v = Gamma(p) - GammaIncompleteComplemented(p, Maths.Pow(x, n));
                return w * Gamma(n) * v;
            }
            return float.NaN;
        }
        /// <summary>
        /// Returns the value of the generalized error function.
        /// </summary>
        /// <param name="x">Value (0, +inf)</param>
        /// <param name="n">Order [0, +inf)</param>
        /// <returns>Value</returns>
        public static Complex32 Gerf(Complex32 x, int n)
        {
            if (n < 0)
                return Complex32.NaN; // singular by your convention

            // E0, E1, E2
            if (n == 0)
                return x / (Maths.E * SQRT_PI);
            if (n == 1)
                return (Complex32.One - Maths.Exp(-x)) / SQRT_PI;
            if (n == 2)
                return Erf(x); // your complex erf

            // General case: En(x) = Γ(n)/√π * γ(1/n, x^n)
            Complex32 xn = Maths.Pow(x, n);                       // x^n, integer power
            Complex32 p = new Complex32(1f / n, 0f);          // 1/n
            Complex32 lower = GammaIncomplete(p, xn);          // γ(1/n, x^n)
            Complex32 gammaN = Gamma(new Complex32(n, 0f));    // Γ(n) = (n-1)!
            return gammaN * lower / SQRT_PI;
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
            return Special.LogGamma(n + 1.0f);
        }
        /// <summary>
        /// Returns the natural logarithm of the factorial of a number log(n!).
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 LogFactorial(Complex32 z)
        {
            return Special.LogGamma(z + Complex32.One);
        }
        /// <summary>
        /// Returns the factorial of a number.
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Value</returns>
        public static double Factorial(float n)
        {
            // check it:
            if (n >= 0 && n == Maths.Round(n))
            {
                // get it from memory
                if (n <= 170)
                {
                    return Special.A000142[(int)n];
                }
                return float.NaN;
            }
            return Special.Gamma(n + 1);
        }
        /// <summary>
        /// Returns the factorial of a number.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Factorial(Complex32 z)
        {
            return Special.Gamma(z + Complex32.One);
        }
        /// <summary>
        /// Returns the decreasing factorial of a number.
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static double FactorialDown(float n, float k)
        {
            return Special.Factorial(n) / Special.Factorial(n - k);
        }
        /// <summary>
        /// Returns the decreasing factorial of a number.
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static Complex32 FactorialDown(Complex32 z, Complex32 k)
        {
            return Special.Gamma(z + Complex32.One) / Special.Gamma(z + Complex32.One - k);
        }
        /// <summary>
        /// Returns the increasing factorial of a number (Pohhammer symbol).
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static float FactorialUp(float n, float k)
        {
            if (n == 0)
            {
                return 1.0f;
            }
            return Special.Gamma(n + k) / Special.Gamma(n);
        }
        /// <summary>
        /// Returns the increasing factorial of a number (Pohhammer symbol).
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static Complex32 FactorialUp(Complex32 z, Complex32 k)
        {
            return Special.Gamma(z + k) / Special.Gamma(z);
        }
        #endregion

        #region Binomial function
        /// <summary>
        /// Returns the value of binomial coefficients: C(n, k) = n! / k! / (n-k)! для k > 0.
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static double Binomial(float n, float k)
        {
            if (k < 0)
            {
                return 0;
            }
            return Special.Factorial(n) / Special.Factorial(k) / Special.Factorial(n - k);
        }
        /// <summary>
        /// Returns the value of binomial coefficients: C(n, k) = n! / k! / (n-k)! для k > 0.
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Binomial(Complex32 n, Complex32 k)
        {
            return Special.Factorial(n) / Special.Factorial(k) / Special.Factorial(n - k);
        }
        /// <summary>
        /// Returns the natural logarithm of binomial coefficients: log(C(n, k)) = log(n!) - log(k!) - log(n-k!).
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static float LogBinomial(float n, float k)
        {
            return Special.LogFactorial(n) - Special.LogFactorial(k) - Special.LogFactorial(n - k);
        }
        /// <summary>
        /// Returns the natural logarithm of binomial coefficients: log(C(n, k)) = log(n!) - log(k!) - log(n-k!).
        /// </summary>
        /// <param name="n">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        public static Complex32 LogBinomial(Complex32 n, Complex32 k)
        {
            return Special.LogFactorial(n) - Special.LogFactorial(k) - Special.LogFactorial(n - k);
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
            if (x == 0f) return 0f;

            if (!inverse)
            {
                // erf(x)
                return (float)ErfApprox(x);
            }
            else
            {
                // inverse error function: y = erf^{-1}(x)
                double y = x;

                if (double.IsNaN(y)) return float.NaN;
                if (y >= 1.0) return float.PositiveInfinity;
                if (y <= -1.0) return float.NegativeInfinity;
                if (y == 0.0) return 0f;

                const double a = 0.147;
                double ya = Math.Abs(y);
                double u = Math.Log(1.0 - ya * ya);
                double s = 2.0 / (Math.PI * a) + 0.5 * u;
                double w = Maths.CopySign((float)Math.Sqrt(Math.Max(0.0f, Math.Sqrt(s * s - u / a) - s)), (float)y);

                for (int it = 0; it < 2; it++)
                {
                    double f = ErfApprox(w) - y;
                    double fp = 2.0 / Math.Sqrt(Math.PI) * Math.Exp(-w * w);
                    double fpp = -2.0 * w * fp;
                    double denom = 2.0 * fp * fp - f * fpp;
                    w -= (denom == 0.0) ? (f / fp) : (2.0 * f * fp / denom);
                }

                return (float)w;
            }
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="inverse">Reverse function or not</param>
        /// <returns>Value</returns>
        public static Complex32 Erf(Complex32 x, bool inverse)
        {
            if (x.Real == 0f && x.Imag == 0f) return Complex32.Zero;

            Complex z = new Complex(x.Real, x.Imag);

            if (!inverse)
            {
                // erf(z)
                return ErfApprox(z);
            }
            else
            {
                // inverse error function: y = erf^{-1}(x)
                if (x.Imag == 0f) return new Complex32(Erf(x.Real, true), 0f);

                const double a = 0.147;
                Complex u = Complex.Log(Complex.One - z * z);
                Complex s = 2.0 / (Math.PI * a) + 0.5 * u;
                Complex rad = Complex.Sqrt(s * s - u / a) - s;
                Complex W = Complex.Sqrt(rad);

                if (W == Complex.Zero)
                {
                    W = z * (Math.Sqrt(Math.PI) / 2.0);
                }
                else
                {
                    if (z.Real < 0) W = -W;
                }

                // f(W)=erf(W)-Z
                for (int it = 0; it < 3; it++)
                {
                    Complex f = ErfApprox(W) - z;
                    Complex fp = 2.0 / Math.Sqrt(Math.PI) * Complex.Exp(-(W * W));      // f'
                    Complex fpp = -2.0 * W * fp;                                        // f''

                    Complex denom = 2.0 * fp * fp - f * fpp;
                    double den2 = denom.Real * denom.Real + denom.Imaginary * denom.Imaginary;
                    Complex step = (den2 < 1e-30) ? (f / fp) : (2.0 * f * fp / denom);
                    W -= step;

                    if (Complex.Abs(step) <= 1e-16 * (1.0 + Complex.Abs(W))) break;
                }

                return new Complex32((float)W.Real, (float)W.Imaginary);
            }
        }
        /// <summary>
        /// Returns the value of the imaginary error function.
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static float Erfi(float x)
        {
            if (x == 0f) return 0f;

            const double DBL_OVER_X = 26.658987628;
            double ax = Math.Abs((double)x);

            if (ax <= 6.0)
            {
                var iz = new Complex(0.0, x);
                var w = ErfApprox(iz);                // erf(i x)
                double val = w.Imaginary;             // -i*(a+ib) = b - ia → Re = b
                if (Math.Abs(val) > float.MaxValue) return Math.Sign(x) * float.PositiveInfinity;
                return (float)val;
            }

            if (ax > DBL_OVER_X) return x > 0 ? float.PositiveInfinity : float.NegativeInfinity;

            // erfi(x) ≈ e^{x^2} / (√π x) * ( 1 + 1/(2x^2) + 3/(4x^4) + 15/(8x^6) + 105/(16x^8) )
            double invx = 1.0 / ax;
            double invx2 = invx * invx;
            double S = 1.0, c = 1.0;

            for (int n = 0; n < 6; n++)
            {
                c *= (2.0 * n + 1.0) * 0.5 * invx2;
                S += c;
                if (Math.Abs(c) < 1e-18 * Math.Abs(S)) break;
            }

            double valp = Math.Exp(ax * ax) * (S / (Math.Sqrt(Math.PI) * ax));
            double valr = x > 0 ? valp : -valp;

            return (float)valr;
        }
        /// <summary>
        /// Returns the value of the imaginary error function.
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static Complex32 Erfi(Complex32 x)
        {
            if (x.Real == 0f && x.Imag == 0f) return Complex32.Zero;

            // erfi(z) = -i * erf(i z)
            var Z = new Complex(x.Real, x.Imag);
            var iZ = new Complex(-Z.Imaginary, Z.Real);  // i*z
            var w = ErfApprox(iZ);                       // erf(i z)

            // -i * (a + i b) = b - i a
            var erfi = new Complex(w.Imaginary, -w.Real);
            return new Complex32((float)erfi.Real, (float)erfi.Imaginary);
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
            return Erf((x - a) / b);
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
            return Erf((x - a) / b);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static float Erfc(float x)
        {
            return 1.0f - Erf(x);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Value</returns>
        public static Complex32 Erfc(Complex32 x)
        {
            return 1.0f - Erf(x);
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
            return 1.0f - Erf(x, a, b);
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
            return 1.0f - Erf(x, a, b);
        }

        #region Private methods (helpers)

        /// <summary>
        /// Approximates the error function for complex arguments.
        /// </summary>
        /// <param name="z">Complex value</param>
        private static Complex ErfApprox(Complex z)
        {
            double az = Complex.Abs(z);

            // |z|: Abramowitz–Stegun 7.1.26
            if (az <= 6.0)
            {
                const double p = 0.3275911;
                const double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
                Complex t = Complex.One / (Complex.One + p * z);
                Complex poly = ((((a5 * t + a4) * t) + a3) * t + a2) * t + a1;
                return Complex.One - poly * t * Complex.Exp(-(z * z));
            }
            else
            {
                // erf(z) ≈ 1 − e^{−z²}/(√π z) * (1 − 1/(2 z²) + 3/(4 z⁴) − 15/(8 z⁶))
                Complex invz = Complex.One / z;
                Complex invz2 = invz * invz;
                Complex S = Complex.One
                          - 0.5 * invz2
                          + 0.75 * (invz2 * invz2)
                          - 1.875 * (invz2 * invz2 * invz2); // 15/8 = 1.875
                return Complex.One - Complex.Exp(-(z * z)) * (invz / Math.Sqrt(Math.PI)) * S;
            }
        }

        static double ErfApprox(double x)
        {
            // |z|: Abramowitz–Stegun 7.1.26
            double ax = Math.Abs(x);
            const double p = 0.3275911;
            const double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
            double t = 1.0 / (1.0 + p * ax);
            double poly = ((((a5 * t + a4) * t) + a3) * t + a2) * t + a1;
            double erfp = 1.0 - poly * t * Math.Exp(-ax * ax);
            return Maths.CopySign((float)erfp, (float)x);
        }

        #endregion

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
            if (positive)
            {
                // D+ function:
                float p = Special.SQRT_PI / 2.0f;
                float v = Maths.Exp(-x * x);
                float erfi = Special.Erfi(x);
                return p * v * erfi;
            }
            // D- function:
            float y = x * x;
            float g = SQRT_PI / 2.0f;
            float e = Special.Erf(x);
            float d = Maths.Exp(y);
            return g * d * e;
        }
        /// <summary>
        /// Returns the value of the D- / D + Dawson function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="positive">D- or D+</param>
        /// <returns>Value</returns>
        public static Complex32 Dawson(Complex32 x, bool positive)
        {
            if (positive)
            {
                // D+ function:
                Complex32 p = Special.SQRT_PI / 2.0f;
                Complex32 v = Maths.Exp(-x * x);
                Complex32 erfi = Special.Erfi(x);
                return p * v * erfi;
            }
            // D- function:
            Complex32 y = x * x;
            Complex32 g = SQRT_PI / 2.0f;
            Complex32 e = Special.Erf(x);
            Complex32 d = Maths.Exp(y);
            return g * d * e;
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
            float re = Maths.Exp(-x * x);
            float im = (float)(2.0 / Math.Sqrt(Math.PI)) * Special.Dawson(x, true);
            return new Complex32(re, im);
        }
        /// <summary>
        /// Returns the value of the Faddeeva function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Faddeeva(Complex32 x)
        {
            return Maths.Exp(-x * x) * Special.Erfc(-x * Complex32.I);
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
            if (inverse)
            {
                return Maths.Sqrt(2) * Special.Erf(1 - 2 * x, true);
            }
            return 0.5f * Special.Erfc(x / Maths.Sqrt2);
        }
        /// <summary>
        /// Returns the value of a Q function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="inverse">Inverse function or not</param>
        /// <returns>Value</returns>
        public static Complex32 Q(Complex32 x, bool inverse = false)
        {
            if (inverse)
            {
                return Maths.Sqrt(2) * Special.Erf(1 - 2 * x, true);
            }
            return 0.5f * Special.Erfc(x / Maths.Sqrt2);
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
            // for all z = 0:
            if (z == 0)
                return 1;

            // Properties:
            double s = 1.0f;
            double m = 1.0f;
            double pa = 1, pb = 1, pc = 1;
            double t;
            float eps = 1e-16f;
            int i, j, iterations = 240;

            // Taylor series:
            for (i = 1; i < iterations; i++)
            {
                // Pochhammer symbols:
                j = i - 1;
                pa *= a + j;
                pb *= b + j;
                pc *= c + j;

                // value:
                m *= z / i;
                t = pa * pb * m / pc;

                // stop point:
                if (Math.Abs(t) < eps)
                { break; }
                else { s += t; }
            }

            // result:
            return (float)s;
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
            // z = 0 => 1
            if (z.Real == 0f && z.Imag == 0f) return Complex32.One;

            Complex s = Complex.One;     // partial sum
            Complex m = Complex.One;     // z^n / n!
            Complex pa = Complex.One;    // (a)_n
            Complex pb = Complex.One;    // (b)_n
            Complex pc = Complex.One;    // (c)_n

            float eps = 1e-16f;
            int maxIter = 240;

            for (int i = 1; i < maxIter; i++)
            {
                double jf = i - 1;
                var j = new Complex(jf, 0f);

                pa *= (Complex)a + j;
                pb *= (Complex)b + j;
                pc *= (Complex)c + j;

                m *= (Complex)z / i; // z^i / i!
                Complex t = pa * pb * m / pc;

                if (Maths.Abs(t) < eps) break;
                s += t;
            }
            return s;
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
            // z = 0 => 1
            if (z == 0f) return 1f;

            bool aNaN = float.IsNaN(a);
            bool bNaN = float.IsNaN(b);

            // 0F0(;;z) = exp(z)  and  1F1(a;a;z) = exp(z)
            if ((aNaN && bNaN) || a == b)
                return Maths.Exp(z);

            // 1F0(a;;z) = (1 - z)^(-a)  (principal branch)
            if (bNaN)
                return Maths.Exp(-a * Maths.Log(1f - z));

            // pole guard: b ∈ {0, -1, -2, ...} (real axis)
            {
                float k = Maths.Round(b);
                if (b <= 0f && Maths.Abs(b - k) < 1e-6f)
                    return float.NaN;
            }

            const float eps = 1e-16f;     // relative stop; adjust to 1e-8..1e-12
            const int iterations = 240;
            const float tiny = 1e-30f;    // protect from near-zero denominator

            double s = 1f;  // partial sum
            double t = 1f;  // current term

            // 0F1(;b;z): t_{n+1} = t_n * z / ((b+n)(n+1))
            if (aNaN)
            {
                for (int n = 0; n < iterations; n++)
                {
                    double bn = b + n;
                    double denom = bn * (n + 1f);

                    if (Maths.Abs(denom) < tiny) return float.NaN; // pole/underflow

                    t *= z / denom;
                    s += t;

                    if (Maths.Abs(t) < eps * (1f + Maths.Abs(s))) break;
                }
                return (float)s;
            }

            // 1F1(a;b;z): t_{n+1} = t_n * (a+n)/(b+n) * z/(n+1)
            // Kummer transform helps convergence if z < 0
            if (z < 0f)
                return Maths.Exp(z) * Hypergeom(b - a, b, -z);

            for (int n = 0; n < iterations; n++)
            {
                double an = a + n;
                double bn = b + n;
                double denom = (n + 1f) * bn;

                if (Maths.Abs(denom) < tiny) return float.NaN; // pole/underflow

                t *= an * z / denom;
                s += t;

                if (Maths.Abs(t) < eps * (1f + Maths.Abs(s))) break;
            }
            return (float)s;
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
            // z = 0 => 1
            if (z.Real == 0f && z.Imag == 0f) return Complex32.One;

            bool aNaN = Complex32.IsNaN(a);
            bool bNaN = Complex32.IsNaN(b);

            // 0F0(;;z) = exp(z)
            if (aNaN && bNaN || (a.Real == b.Real && a.Imag == b.Imag))
                return Maths.Exp(z);

            // 1F0(a;;z) = (1 - z)^(-a)  (principal branch)
            if (bNaN)
                return Maths.Exp(-a * Maths.Log(Complex32.One - z));

            // pole guard: b ∈ {0, -1, -2, ...} (real axis)
            if (b.Imag == 0f)
            {
                float br = b.Real, k = Maths.Round(br);
                if (br <= 0f && Maths.Abs(br - k) < 1e-6f)
                    return Complex32.NaN;
            }

            const float eps = 1e-16f;
            const int iterations = 240;
            const float tiny = 1e-30f;

            Complex s = Complex.One;   // partial sum
            Complex t = Complex.One;   // current term

            // 0F1(;b;z): t_{n+1} = t_n * z / ((b+n)(n+1))
            if (aNaN)
            {
                for (int n = 0; n < iterations; n++)
                {
                    Complex bn = (Complex)b + new Complex(n, 0f);
                    Complex denom = bn * new Complex(n + 1f, 0f);

                    if (Maths.Abs(denom) < tiny) return Complex32.NaN; // pole/underflow

                    t *= (Complex)z / denom;
                    s += t;

                    if (Maths.Abs(t) < eps * (1f + Maths.Abs(s))) break;
                }
                return s;
            }

            // 1F1(a;b;z): t_{n+1} = t_n * (a+n)/(b+n) * z/(n+1)
            // Kummer transform helps convergence if Re(z) < 0
            if (z.Real < 0f)
                return Maths.Exp(z) * Hypergeom(b - a, b, -z);

            for (int n = 0; n < iterations; n++)
            {
                Complex an = (Complex)a + new Complex(n, 0f);
                Complex bn = (Complex)b + new Complex(n, 0f);
                Complex denom = new Complex(n + 1f, 0f) * bn;

                if (Maths.Abs(denom) < tiny) return Complex32.NaN; // pole/underflow

                t *= an * (Complex)z / denom;
                s += t;

                if (Maths.Abs(t) < eps * (1f + Maths.Abs(s))) break;
            }
            return s;
        }
        #endregion

        #region Beta functions
        /// <summary>
        /// Returns the value of the beta function: B(a, b) = Г(a) * Г(b) / Г(ab).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static float Beta(float a, float b)
        {
            return Special.Gamma(a) * Special.Gamma(b) / Special.Gamma(a + b);
        }
        /// <summary>
        /// Returns the value of the beta function: B(a, b) = Г(a) * Г(b) / Г(ab).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Beta(Complex32 a, Complex32 b)
        {
            return Special.Gamma(a) * Special.Gamma(b) / Special.Gamma(a + b);
        }
        /// <summary>
        /// Returns the value of the beta function: B(m, n) = (m - 1)! * (n - 1)! / (m + n - 1)!.
        /// </summary>
        /// <param name="m">Integer number</param>
        /// <param name="n">Integer number</param>
        /// <returns>Value</returns>
        public static double Beta(int m, int n)
        {
            return Special.Factorial(m - 1) * Special.Factorial(n - 1) / Special.Factorial(m + n - 1);
        }
        /// <summary>
        /// Returns the value of a derivative beta function: B'(a, b).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static float BetaDerivative(float a, float b)
        {
            return Special.Beta(a, b) * (Special.DiGamma(a) - Special.DiGamma(a + b));
        }
        /// <summary>
        /// Returns the value of a derivative beta function: B'(a, b).
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static Complex32 BetaDerivative(Complex32 a, Complex32 b)
        {
            return Special.Beta(a, b) * (Special.DiGamma(a) - Special.DiGamma(a + b));
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
            if (float.IsNaN(a) || float.IsNaN(b) || float.IsNaN(x)) return float.NaN;

            // Domain commonly assumed for the real-valued branch:
            if (x <= 0f) return 0f;
            if (x >= 1f)
            {
                // Complete beta B(a,b)
                return Special.Gamma(a) * Special.Gamma(b) / Special.Gamma(a + b);
            }

            // Poles of Beta when a or b is a non-positive integer
            if ((a <= 0f && a == Maths.Round(a)) || (b <= 0f && b == Maths.Round(b)))
                return float.NaN;

            // Optional symmetry to improve conditioning near x≈1:
            // B_x(a,b) = B(a,b) - B_{1-x}(b,a)
            if (x > 0.5f)
            {
                float B = Special.Gamma(a) * Special.Gamma(b) / Special.Gamma(a + b);
                return B - BetaIncomplete(b, a, 1f - x);
            }

            // Hypergeometric representation
            float pref = Maths.Exp(a * Maths.Log(x)) / a;          // x^a / a
            float hyp = Special.Hypergeom(a, 1f - b, a + 1f, x);  // 2F1(a,1-b; a+1; x)
            return pref * hyp;
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
            // B_0(a,b) = 0
            if (x.Real == 0f && x.Imag == 0f) return Complex32.Zero;

            // B_1(a,b) = B(a,b) = Γ(a)Γ(b)/Γ(a+b)
            if (x.Imag == 0f && x.Real == 1f)
                return Special.Gamma(a) * Special.Gamma(b) / Special.Gamma(a + b);

            // Principal-branch power and Gauss 2F1
            Complex32 pref = Maths.Exp(a * Maths.Log(x)) / a;                   // x^a / a
            Complex32 hyp = Special.Hypergeom(a, Complex32.One - b, a + Complex32.One, x); // 2F1(a,1-b; a+1; x)
            return pref * hyp;
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
            return Special.BetaIncomplete(a, b, x) / Special.Beta(a, b);
        }
        /// <summary>
        /// Returns the value of a log-beta function.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static float LogBeta(float a, float b)
        {
            return Special.LogGamma(a) + Special.LogGamma(b) - Special.LogGamma(a + b);
        }
        /// <summary>
        /// Returns the value of a log-beta function.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static Complex32 LogBeta(Complex32 a, Complex32 b)
        {
            return Special.LogGamma(a) + Special.LogGamma(b) - Special.LogGamma(a + b);
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
            if (x == 0f) return float.NegativeInfinity;
            if (x < 0f) return float.NaN;

            double xd = x;
            const double XSW = 16.0;

            if (xd <= XSW)
            {
                // Ci(x) = γ + ln x + ∫_0^1 (cos(x t) - 1)/t dt
                double integral = GaussLegendre01(t => (Math.Cos(xd * t) - 1.0) / t);
                return (float)(Maths.Gamma + Math.Log(xd) + integral);
            }
            else
            {
                // S0 = Σ (-1)^n (2n)! / x^{2n+1}, S1 = Σ (-1)^n (2n+1)! / x^{2n+2}
                double inv = 1.0 / xd, inv2 = inv * inv;
                double s0 = inv * (1.0 - 2.0 * inv2 + 24.0 * inv2 * inv2 - 720.0 * inv2 * inv2 * inv2);
                double s1 = inv * inv * (1.0 - 6.0 * inv2 + 120.0 * inv2 * inv2 - 5040.0 * inv2 * inv2 * inv2);
                double sinx = Math.Sin(xd), cosx = Math.Cos(xd);
                return (float)(sinx * s0 - cosx * s1);
            }

            static double GaussLegendre01(Func<double, double> f, int panels = 1)
            {
                double[] u = new double[] {
                    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                    -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                     0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                     0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
                };

                double[] w = new double[] {
                    0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                    0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                    0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                    0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
                };

                double sum = 0.0, h = 1.0 / panels;
                for (int p = 0; p < panels; p++)
                {
                    double a = p * h, b = (p + 1) * h;
                    double half = 0.5 * (b - a), mid = 0.5 * (a + b);
                    double loc = 0.0;
                    for (int i = 0; i < 16; i++)
                    {
                        double t = mid + half * u[i];
                        loc += w[i] * f(t);
                    }
                    sum += half * loc;
                }
                return sum;
            }
        }
        /// <summary>
        /// Returns the value of the integral cosine.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Ci(Complex32 z)
        {
            if (z.Real == 0f && z.Imag == 0f)
                return new Complex32(float.NegativeInfinity, 0f);

            var Z = new Complex(z.Real, z.Imag);
            double r = Complex.Abs(Z);
            const double XSW = 16.0;

            if (r <= XSW)
            {
                // Ci(z) = γ + Log(z) + ∫_0^1 (cos(z t) - 1)/t dt  (principal log)
                var integral = GaussLegendre01C(t => Complex.Cos(Z * t) - Complex.One);
                var val = (Maths.Gamma + Complex.Log(Z)) + integral;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
            else
            {
                var inv = Complex.One / Z;
                var inv2 = inv * inv;
                var s0 = inv * (1.0 - 2.0 * inv2 + 24.0 * inv2 * inv2 - 720.0 * inv2 * inv2 * inv2);
                var s1 = inv2 * (1.0 - 6.0 * inv2 + 120.0 * inv2 * inv2 - 5040.0 * inv2 * inv2 * inv2);
                var sinz = Complex.Sin(Z);
                var cosz = Complex.Cos(Z);
                var val = sinz * s0 - cosz * s1;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }

            static Complex GaussLegendre01C(Func<double, Complex> f, int panels = 1)
            {
                double[] u = new double[] {
                    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                    -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                     0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                     0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
                };
                double[] w = new double[] {
                    0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                    0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                    0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                    0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
                };
                Complex sum = Complex.Zero;
                double h = 1.0 / panels;
                for (int p = 0; p < panels; p++)
                {
                    double a = p * h, b = (p + 1) * h;
                    double half = 0.5 * (b - a), mid = 0.5 * (a + b);
                    Complex loc = Complex.Zero;
                    for (int i = 0; i < 16; i++)
                    {
                        double t = mid + half * u[i];
                        loc += w[i] * (f(t) / t);
                    }
                    sum += half * loc;
                }
                return sum;
            }
        }
        /// <summary>
        /// Returns the value of the integral sine.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Si(float x)
        {
            if (x < 0f) return -Special.Si(-x);

            double xd = x;
            const double XSW = 16.0;

            if (xd <= XSW)
            {
                // Si(x) = ∫_0^1 sin(x t)/t dt
                double integral = Ci_GaussLegendre01(t => Math.Sin(xd * t) / t);
                return (float)integral;
            }
            else
            {
                // Si ≈ π/2 − cos x * S0 − sin x * S1
                double inv = 1.0 / xd, inv2 = inv * inv;
                double s0 = inv * (1.0 - 2.0 * inv2 + 24.0 * inv2 * inv2 - 720.0 * inv2 * inv2 * inv2);
                double s1 = inv * inv * (1.0 - 6.0 * inv2 + 120.0 * inv2 * inv2 - 5040.0 * inv2 * inv2 * inv2);
                return (float)(Math.PI * 0.5 - Math.Cos(xd) * s0 - Math.Sin(xd) * s1);
            }

            static double Ci_GaussLegendre01(Func<double, double> f, int panels = 1)
            {
                double[] u = new double[] {
                    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                    -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                     0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                     0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
                };

                double[] w = new double[] {
                    0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                    0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                    0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                    0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
                };

                double sum = 0.0, h = 1.0 / panels;
                for (int p = 0; p < panels; p++)
                {
                    double a = p * h, b = (p + 1) * h;
                    double half = 0.5 * (b - a), mid = 0.5 * (a + b);
                    double loc = 0.0;
                    for (int i = 0; i < 16; i++)
                    {
                        double t = mid + half * u[i];
                        loc += w[i] * f(t);
                    }
                    sum += half * loc;
                }
                return sum;
            }
        }
        /// <summary>
        /// Returns the value of the integral sine.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Si(Complex32 z)
        {
            var Z = new Complex(z.Real, z.Imag);
            double r = Complex.Abs(Z);
            const double XSW = 16.0;

            if (r <= XSW)
            {
                // Si(z) = ∫_0^1 sin(z t)/t dt
                var integral = GL01C(t => Complex.Sin(Z * t) / t);
                return new Complex32((float)integral.Real, (float)integral.Imaginary);
            }
            else
            {
                var inv = Complex.One / Z;
                var inv2 = inv * inv;
                var s0 = inv * (1.0 - 2.0 * inv2 + 24.0 * inv2 * inv2 - 720.0 * inv2 * inv2 * inv2);
                var s1 = inv2 * (1.0 - 6.0 * inv2 + 120.0 * inv2 * inv2 - 5040.0 * inv2 * inv2 * inv2);
                var val = Math.PI * 0.5 - Complex.Cos(Z) * s0 - Complex.Sin(Z) * s1;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }

            static Complex GL01C(Func<double, Complex> f, int panels = 1)
            {
                double[] u = new double[] {
                    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                    -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                     0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                     0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
                };

                double[] w = new double[] {
                    0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                    0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                    0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                    0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
                };

                Complex sum = Complex.Zero;
                double h = 1.0 / panels;
                for (int p = 0; p < panels; p++)
                {
                    double a = p * h, b = (p + 1) * h;
                    double half = 0.5 * (b - a), mid = 0.5 * (a + b);
                    Complex loc = Complex.Zero;
                    for (int i = 0; i < 16; i++)
                    {
                        double t = mid + half * u[i];
                        loc += w[i] * f(t);
                    }
                    sum += half * loc;
                }
                return sum;
            }
        }
        /// <summary>
        /// Returns the value of an integral exponential function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Ei(float x)
        {
            if (x == 0f) return float.NegativeInfinity;

            double xd = x;
            const double XASY = 25.0;

            if (Math.Abs(xd) >= XASY)
            {
                double inv = 1.0 / xd, inv2 = inv * inv;
                double poly = 1.0 + inv + 2.0 * inv2 + 6.0 * inv2 * inv + 24.0 * inv2 * inv2;
                double val = Math.Exp(xd) * inv * poly;

                return (float)val;
            }

            // Ei(x) = γ + ln|x| + ∫_0^1 (e^{x t} - 1)/t dt  (|arg x|<π)
            double integral = GL01(t => (Math.Exp(xd * t) - 1.0) / t);
            return (float)(Maths.Gamma + Math.Log(Math.Abs(xd)) + integral);

            static double GL01(Func<double, double> f, int panels = 1)
            {
                double[] u = new double[] {
                    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                    -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                     0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                     0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
                };
                double[] w = new double[] {
                    0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                    0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                    0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                    0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
                };

                double sum = 0.0;
                for (int i = 0; i < 16; i++)
                {
                    double t = 0.5 * (u[i] + 1.0);     // [0,1]
                    sum += w[i] * f(t);
                }
                return 0.5 * sum;
            }
        }
        /// <summary>
        /// Returns the value of an integral exponential function.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Ei(Complex32 z)
        {
            if (z.Real == 0f && z.Imag == 0f)
                return new Complex32(float.NegativeInfinity, 0f);

            var Z = new Complex(z.Real, z.Imag);
            double r = Complex.Abs(Z);
            const double XASY = 25.0;

            if (r >= XASY)
            {
                //  Ei(z) ~ e^z / z * (1 + 1/z + 2!/z^2 + 3!/z^3 + 4!/z^4), |arg z|<π
                var inv = Complex.One / Z;
                var inv2 = inv * inv;
                var poly = Complex.One + inv + 2.0 * inv2 + 6.0 * inv2 * inv + 24.0 * inv2 * inv2;
                var val = Complex.Exp(Z) * inv * poly;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
            else
            {
                // Ei(z) = γ + Log(z) + ∫_0^1 (e^{z t} - 1)/t dt
                var integral = GL01C(t => (Complex.Exp(Z * t) - Complex.One) / t);
                var val = (Maths.Gamma + Complex.Log(Z)) + integral;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }

            static Complex GL01C(Func<double, Complex> f)
            {
                double[] u = new double[] {
                    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                    -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                     0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                     0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
                };
                double[] w = new double[] {
                    0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                    0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                    0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                    0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
                };

                Complex sum = Complex.Zero;
                for (int i = 0; i < 16; i++)
                {
                    double t = 0.5 * (u[i] + 1.0); // [0,1]
                    sum += w[i] * f(t);
                }
                return 0.5 * sum;
            }
        }
        /// <summary>
        /// Returns the value of the integral logarithm.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Li(float x)
        {
            if (x < 0f) return float.NaN;
            return Ei(Maths.Log(x));
        }
        /// <summary>
        /// Returns the value of the integral logarithm.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Li(Complex32 z)
        {
            return Ei(Maths.Log(z));
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
            // J_n(0) = δ_{n0}
            if (x == 0f) return (a == 0) ? 1f : 0f;

            int n = a < 0 ? -a : a;
            double signOrder = (a < 0 && ((n & 1) == 1)) ? -1.0 : 1.0;
            const double XSW = 20.0;

            // J_n(-x) = (-1)^n J_n(x)
            double xabs = x;
            double signArg = 1.0;

            if (x < 0.0) { xabs = -x; signArg = ((n & 1) == 1) ? -1.0 : 1.0; }

            if (xabs < XSW)
            {
                double val = GL0Pi(theta => Math.Cos(n * theta - xabs * Math.Sin(theta)));
                return (float)(signOrder * signArg * (val / Math.PI));
            }
            else
            {
                double mu = 4.0 * n * n;
                double phi = xabs - 0.5 * n * Math.PI - Math.PI / 4.0;
                double inv = 1.0 / xabs;

                double c0 = Math.Sqrt(2.0 / (Math.PI * xabs));
                double A1 = (mu - 1.0) * 0.125 * inv;                       // (μ-1)/(8x)
                double A2 = (mu - 1.0) * (mu - 9.0) * (inv * inv) / 128.0;  // ((μ-1)(μ-9))/(128 x^2)

                double val = c0 * (Math.Cos(phi) - A1 * Math.Sin(phi) - A2 * Math.Cos(phi));
                return (float)(signOrder * signArg * val);
            }
        }
        /// <summary>
        /// Returns the value of a Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 J(Complex32 x, int a)
        {
            if (x.Real == 0f && x.Imag == 0f)
                return (a == 0) ? Complex32.One : Complex32.Zero;

            int n = a < 0 ? -a : a;
            double signOrder = (a < 0 && ((n & 1) == 1)) ? -1.0 : 1.0;
            const double XSW = 20.0;

            // J_n(-x)=(-1)^n J_n(x)
            if (x.Imag == 0f)
            {
                return J(x.Real, a);
            }

            Complex Z = new Complex(x.Real, x.Imag);
            double r = Complex.Abs(Z);

            if (r < XSW)
            {
                // J_n(z) = (1/π) ∫_0^π cos(n t − z sin t) dt
                Complex sum = GL0PiC(theta => Complex.Cos(n * theta - Z * Math.Sin(theta)));
                Complex val = sum / Math.PI * signOrder;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
            else
            {
                double mu = 4.0 * n * n;
                double c0 = Math.Sqrt(2.0 / (Math.PI * r));
                double A1 = (mu - 1.0) * 0.125 / r;
                double A2 = (mu - 1.0) * (mu - 9.0) / 128.0 / (r * r);

                Complex phi = Z - (n * Math.PI * 0.5 + Math.PI * 0.25);
                Complex val = c0 * (Complex.Cos(phi) - A1 * Complex.Sin(phi) - A2 * Complex.Cos(phi));
                val *= signOrder;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
        }

        /// <summary>
        /// Returns the value of a Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Y(float x, int a)
        {
            if (x <= 0f) return float.NaN;

            int n = a < 0 ? -a : a;
            double sgn = (a < 0 && ((n & 1) == 1)) ? -1.0 : 1.0;

            double ax = x;
            const double XSW = 16.0;

            if (ax <= XSW)
            {
                // Y_n(x) = (1/π) ∫_0^π sin(x sin t - n t) dt - (1/π) ∫_0^∞ e^{-x sinh u - n u} du
                double Iosc = GL0Pi(t => Math.Sin(x * Math.Sin(t) - n * t)) / Math.PI;
                double T = Math.Max(8.0, Math.Log(40.0 / x + 1.0) + 2.0);
                double Idc = GL0T(u => Math.Exp(-x * Math.Sinh(u) - n * u), T) / Math.PI;
                double val = Iosc - Idc;
                return (float)(sgn * val);
            }
            else
            {
                double mu = 4.0 * n * n;
                double phi = ax - 0.5 * n * Math.PI - Math.PI / 4.0;
                double inv = 1.0 / ax;
                double c0 = Math.Sqrt(2.0 / (Math.PI * ax));
                double A1 = (mu - 1.0) * 0.125 * inv;
                double A2 = (mu - 1.0) * (mu - 9.0) / 256.0 * inv * inv;

                double val = c0 * (Math.Sin(phi) + A1 * Math.Cos(phi) - A2 * Math.Sin(phi));
                return (float)(val * sgn);
            }
        }
        /// <summary>
        /// Returns the value of a Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Y(Complex32 x, int a)
        {
            if (x.Real == 0f && x.Imag == 0f) return Complex32.NaN;

            int n = a < 0 ? -a : a;
            double sgn = (a < 0 && ((n & 1) == 1)) ? -1.0 : 1.0;

            Complex Z = new Complex(x.Real, x.Imag);
            double r = Complex.Abs(Z);
            const double XSW = 16.0;

            if (r <= XSW && x.Real > 0f)
            {
                // Y_n(z) = (1/π)[ ∫_0^π sin(z sin t - n t) dt - ∫_0^∞ e^{-z sinh u - n u} du ]
                Complex Iosc = GL0PiC(t => Complex.Sin(Z * Math.Sin(t) - n * t)) / Math.PI;
                double T = Math.Max(8.0, Math.Log(40.0 / (r + 1e-12) + 1.0) + 2.0);
                Complex Idc = GL0TC(u => Complex.Exp(-Z * Math.Sinh(u) - n * u), T) / Math.PI;
                Complex val = Iosc - Idc;
                val *= sgn;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
            else
            {
                double mu = 4.0 * n * n;
                double ax = r;
                double c0 = Math.Sqrt(2.0 / (Math.PI * ax));
                double A1 = (mu - 1.0) * 0.125 / ax;
                double A2 = (mu - 1.0) * (mu - 9.0) / 256.0 / (ax * ax);
                Complex cosφ = Complex.Cos(Z - (n * Math.PI * 0.5 + Math.PI * 0.25));
                Complex sinφ = Complex.Sin(Z - (n * Math.PI * 0.5 + Math.PI * 0.25));
                Complex val = c0 * (sinφ + A1 * cosφ - A2 * sinφ);
                val *= sgn;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
        }

        /// <summary>
        /// Returns the value of the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float I(float x, int a)
        {
            if (x == 0f) return (a == 0) ? 1f : 0f;

            int n = a < 0 ? -a : a; // I_{-n} = I_n

            double ax = Math.Abs((double)x);
            const double XSW = 20.0;

            if (ax <= XSW)
            {
                // I_n(x) = (1/π) ∫_0^π e^{x cos t} cos(n t) dt
                double val = GL0Pi(t => Math.Exp(x * Math.Cos(t)) * Math.Cos(n * t)) / Math.PI;
                return (float)val;
            }
            else
            {
                // I_n(x) ~ e^x / √(2πx) [1 - (μ-1)/(8x) + ((μ-1)(μ-9))/(2!·128 x^2)]
                double mu = 4.0 * n * n;
                double inv = 1.0 / ax;
                double c0 = Math.Exp(ax) / Math.Sqrt(2.0 * Math.PI * ax);
                double P1 = -(mu - 1.0) * 0.125 * inv;
                double P2 = ((mu - 1.0) * (mu - 9.0)) / 256.0 * inv * inv;
                double val = c0 * (1.0 + P1 + P2);
                return (float)val;
            }
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 I(Complex32 x, int a)
        {
            if (x.Real == 0f && x.Imag == 0f)
                return (a == 0) ? Complex32.One : Complex32.Zero;

            int n = a < 0 ? -a : a; // I_{-n} = I_n
            Complex Z = new Complex(x.Real, x.Imag);
            double r = Complex.Abs(Z);
            const double XSW = 20.0;

            if (r <= XSW)
            {
                // I_n(z) = (1/π) ∫_0^π e^{z cos t} cos(n t) dt
                Complex val = GL0PiC(t => Complex.Exp(Z * Math.Cos(t)) * Math.Cos(n * t)) / Math.PI;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
            else
            {
                double mu = 4.0 * n * n;
                double ax = r;
                double c0 = Math.Exp(ax) / Math.Sqrt(2.0 * Math.PI * ax);
                double P1 = -(mu - 1.0) * 0.125 / ax;
                double P2 = (mu - 1.0) * (mu - 9.0) / 256.0 / (ax * ax);
                Complex scale = Complex.Exp(Z - ax);
                Complex val = scale * (c0 * (1.0 + P1 + P2));
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
        }

        /// <summary>
        /// Returns the value of the modified Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float K(float x, int a)
        {
            if (x <= 0f) return float.NaN;

            int n = a < 0 ? -a : a; // K_{-n} = K_n

            const double XSW = 20.0;

            if (x <= XSW)
            {
                // K_n(x) = ∫_0^∞ e^{-x cosh t} cosh(n t) dt
                double T = Math.Max(8.0, Math.Log(40.0 / x + 1.0) + 2.0);
                double val = GL0T(t => Math.Exp(-x * Math.Cosh(t)) * Math.Cosh(n * t), T);
                return (float)val;
            }
            else
            {
                // K_n(x) ~ √(π/(2x)) e^{-x} [1 + (μ-1)/(8x) + ((μ-1)(μ-9))/(2!·128 x^2)]
                double mu = 4.0 * n * n;
                double inv = 1.0 / x;
                double c0 = Math.Sqrt(Math.PI / (2.0 * x)) * Math.Exp(-x);
                double P1 = (mu - 1.0) * 0.125 * inv;
                double P2 = ((mu - 1.0) * (mu - 9.0)) / 256.0 * inv * inv;
                double val = c0 * (1.0 + P1 + P2);
                return (float)val;
            }
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 K(Complex32 x, int a)
        {
            if (x.Real == 0f && x.Imag == 0f) return Complex32.NaN;

            int n = a < 0 ? -a : a; // K_{-n} = K_n
            Complex Z = new Complex(x.Real, x.Imag);
            double r = Complex.Abs(Z);
            const double XSW = 20.0;

            if (x.Real > 0f && r <= XSW)
            {
                double T = Math.Max(8.0, Math.Log(40.0 / (r + 1e-12) + 1.0) + 2.0);
                Complex val = GL0TC(t => Complex.Exp(-Z * Math.Cosh(t)) * Math.Cosh(n * t), T);
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
            else
            {
                double mu = 4.0 * n * n;
                double inv = 1.0 / r;
                double c0 = Math.Sqrt(Math.PI / (2.0 * r));

                // K_n(z) ~ √(π/(2z)) e^{-z}(1 + …)
                Complex scale = Complex.Exp(-Z);
                double P1 = (mu - 1.0) * 0.125 * inv;
                double P2 = ((mu - 1.0) * (mu - 9.0)) / 256.0 * inv * inv;
                Complex val = scale * (c0 * (1.0 + P1 + P2));
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
        }

        #region Private methods (helpers)

        /// <summary>
        /// Integrates a real function over [0, π] using 16-point Gauss–Legendre quadrature.
        /// </summary>
        /// <param name="f">Function to integrate</param>
        private static double GL0Pi(Func<double, double> f)
        {
            double[] u = new double[] {
                -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                 0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                 0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
            };
            double[] w = new double[] {
                0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
            };
            double half = Math.PI * 0.5, sum = 0.0;
            for (int i = 0; i < 16; i++)
            {
                double t = half * (u[i] + 1.0);
                sum += w[i] * f(t);
            }
            return half * sum;
        }

        /// <summary>
        /// Integrates a complex function over [0, π] using 16-point Gauss–Legendre quadrature.
        /// </summary>
        /// <param name="f">Function to integrate</param>
        private static Complex GL0PiC(Func<double, Complex> f)
        {
            double[] u = new double[] {
                -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                 0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                 0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
            };
            double[] w = new double[] {
                0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
            };
            double half = Math.PI * 0.5;
            Complex sum = Complex.Zero;
            for (int i = 0; i < 16; i++)
            {
                double t = half * (u[i] + 1.0);
                sum += w[i] * f(t);
            }
            return half * sum;
        }

        /// <summary>
        /// Integrates a real function over [0, T] using 16-point Gauss–Legendre quadrature.
        /// </summary>
        /// <param name="f">Function to integrate</param>
        /// <param name="T">Upper integration limit</param>
        private static double GL0T(Func<double, double> f, double T)
        {
            double[] u = new double[] {
                -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                 0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                 0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
            };
            double[] w = new double[] {
                0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
            };
            double half = 0.5 * T, mid = half;
            double sum = 0.0;
            for (int i = 0; i < 16; i++)
            {
                double t = mid + half * u[i];
                sum += w[i] * f(t);
            }
            return half * sum;
        }

        static Complex GL0TC(Func<double, Complex> f, double T)
        {
            double[] u = new double[] {
                -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
                -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744,
                 0.09501250983763744,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
                 0.7554044083550030,   0.8656312023878318,  0.9445750230732326,  0.9894009349916499
            };
            double[] w = new double[] {
                0.027152459411754095, 0.06225352393864789, 0.09515851168249278, 0.12462897125553387,
                0.14959598881657673,  0.16915651939500254, 0.18260341504492359, 0.18945061045506850,
                0.18945061045506850,  0.18260341504492359, 0.16915651939500254, 0.14959598881657673,
                0.12462897125553387,  0.09515851168249278, 0.06225352393864789, 0.027152459411754095
            };
            double half = 0.5 * T, mid = half;
            Complex sum = Complex.Zero;
            for (int i = 0; i < 16; i++)
            {
                double t = mid + half * u[i];
                sum += w[i] * f(t);
            }
            return half * sum;
        }
        #endregion

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
            if (a < 0) return 0f;
            if (x == 0f) return 0f;

            double ax = Math.Abs((double)x);
            double signArg = (x < 0) ? ((a & 1) == 0) ? -1.0 : 1.0 : 1.0;

            const double XSW = 25.0;
            const double BW = 1.0;

            double invNorm0 = 1.0 / (Math.Sqrt(Math.PI) * Special.Gamma(a + 0.5f)); // 1/(√π Γ(a+1/2))
            double invNorm1 = 1.0 / (Math.Sqrt(Math.PI) * Special.Gamma(a - 0.5f)); // 1/(√π Γ(a-1/2))

            double smallVal()
            {
                double pref = 2.0 * Math.Pow(0.5 * ax, a) * invNorm0;
                int panels = Math.Max(2, (int)Math.Ceiling(ax / 4.0));

                double[] u = {
                    -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                     0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363
                };
                double[] w = {
                    0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                    0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763
                };

                double sum = 0.0;
                for (int p = 0; p < panels; p++)
                {
                    double aτ = p * (Math.PI * 0.5) / panels;
                    double bτ = (p + 1) * (Math.PI * 0.5) / panels;
                    double half = 0.5 * (bτ - aτ);
                    double mid = 0.5 * (aτ + bτ);

                    for (int i = 0; i < 8; i++)
                    {
                        double tau = mid + half * u[i];
                        double wt = half * w[i];
                        double s2a = Math.Pow(Math.Sin(tau), 2 * a);
                        double f = Math.Sin(ax * Math.Cos(tau));
                        sum += wt * s2a * f;
                    }
                }
                return pref * sum;
            }

            double largeVal()
            {
                double ya = Special.Y((float)ax, a); // Y_a(|x|)
                double lead0 = Math.Pow(0.5 * ax, a - 1) * invNorm0;
                double lead1 = 4.0 / 3.0 * Math.Pow(0.5 * ax, a - 3) * invNorm1;
                return ya + (lead0 - lead1);
            }

            double valSmall = smallVal();
            double valLarge = largeVal();
            double val;

            if (ax <= XSW - BW) val = valSmall;
            else if (ax >= XSW + BW) val = valLarge;
            else
            {
                double t = (ax - (XSW - BW)) / (2.0 * BW);
                double wBlend = 0.5 * (1.0 - Math.Cos(Math.PI * t));
                val = (1.0 - wBlend) * valSmall + wBlend * valLarge;
            }

            return (float)(signArg * val);
        }
        /// <summary>
        /// Returns the value of the Struve function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static Complex32 H(Complex32 x, int a)
        {
            if (a < 0) return Complex32.Zero;
            if (x.Real == 0f && x.Imag == 0f) return Complex32.Zero;
            if (x.Imag == 0f) return new Complex32(H(x.Real, a), 0f);

            var Z = new Complex(x.Real, x.Imag);
            double r = Maths.Abs(x);

            const double XSW = 25.0;
            const double BW = 1.0;

            double invNorm0 = 1.0 / (Math.Sqrt(Math.PI) * Special.Gamma(a + 0.5f)); // 1/(√π Γ(a+1/2))
            double invNorm1 = 1.0 / (Math.Sqrt(Math.PI) * Special.Gamma(a - 0.5f)); // 1/(√π Γ(a-1/2))

            Complex smallVal()
            {
                // pref = 2/(√π Γ) * (Z/2)^a
                var pref = 2.0 * invNorm0 * Complex.Pow(Z * 0.5, a);

                int panels = Math.Max(2, (int)Math.Ceiling(r / 4.0));
                double[] u = {
                    -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                     0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363
                };
                double[] w = {
                    0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                    0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763
                };

                var sum = Complex.Zero;
                for (int p = 0; p < panels; p++)
                {
                    double aτ = (p) * (Math.PI * 0.5) / panels;
                    double bτ = (p + 1) * (Math.PI * 0.5) / panels;
                    double half = 0.5 * (bτ - aτ);
                    double mid = 0.5 * (aτ + bτ);

                    for (int i = 0; i < 8; i++)
                    {
                        double tau = mid + half * u[i];
                        double wt = half * w[i];
                        double s2a = Math.Pow(Math.Sin(tau), 2 * a);
                        var f = Complex.Sin(Z * Math.Cos(tau));
                        sum += wt * s2a * f;
                    }
                }
                return pref * sum;
            }

            Complex largeVal()
            {
                Complex32 ya32 = Special.Y(x, a);
                var ya = new Complex(ya32.Real, ya32.Imag);

                var lead0 = Complex.Pow(Z * 0.5, a - 1) * invNorm0;
                var lead1 = (4.0 / 3.0) * Complex.Pow(Z * 0.5, a - 3) * invNorm1;

                return ya + (lead0 - lead1);
            }

            var valSmall = smallVal();
            var valLarge = largeVal();
            Complex val;

            if (r <= XSW - BW) val = valSmall;
            else if (r >= XSW + BW) val = valLarge;
            else
            {
                double t = (r - (XSW - BW)) / (2.0 * BW);
                double wBlend = 0.5 * (1.0 - Math.Cos(Math.PI * t));
                val = (1.0 - wBlend) * valSmall + wBlend * valLarge;
            }

            return new Complex32((float)val.Real, (float)val.Imaginary);
        }
        /// <summary>
        /// Returns the value of the modified Struve function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="v">Value</param>
        /// <returns>Value</returns>
        public static float L(float x, int v)
        {
            if (v < 0) return 0f;
            if (x == 0f) return 0f;

            double ax = Math.Abs((double)x);
            double signArg = (x < 0) ? (((v & 1) == 0) ? -1.0 : 1.0) : 1.0; // (-1)^{v+1}

            const double XSW = 12.0;
            double gv = Special.Gamma(v + 0.5f);
            double invNorm = 1.0 / (Math.Sqrt(Math.PI) * gv);

            if (ax <= XSW)
            {
                double pref = 2.0 * Math.Pow(0.5 * ax, v) * invNorm;

                double sum = 0.0, half = Math.PI * 0.25;
                double[] u = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                        0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363 };
                double[] w = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                       0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };

                for (int i = 0; i < 8; i++)
                {
                    double tau = half * (u[i] + 1.0);
                    double wt = half * w[i];
                    double s2v = Math.Pow(Math.Sin(tau), 2 * v);
                    double f = Math.Sinh(ax * Math.Cos(tau));   // <--- ax
                    sum += wt * s2v * f;
                }

                double val = pref * sum;
                return (float)(signArg * val);
            }
            else
            {
                double lead = Math.Pow(0.5 * ax, v - 1) * invNorm;
                double iv = Special.I((float)ax, v);   // I(|x|,v); I_n(-x) = (-1)^n I_n(x)
                double val = iv + lead;
                return (float)(signArg * val);
            }
        }
        /// <summary>
        /// Returns the value of the modified Struve function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="v">Value</param>
        /// <returns>Value</returns>
        public static Complex32 L(Complex32 x, int v)
        {
            if (v < 0) return Complex32.Zero;
            if (x.Real == 0f && x.Imag == 0f) return Complex32.Zero;

            if (x.Imag == 0f)
                return new Complex32(L(x.Real, v), 0f);

            var Z = new Complex(x.Real, x.Imag);
            double r = Maths.Abs(x);
            const double XSW = 12.0;

            double gv = Special.Gamma(v + 0.5f);
            double invNorm = 1.0 / (Math.Sqrt(Math.PI) * gv);

            if (r <= XSW)
            {
                var sum = Complex.Zero;
                double half = Math.PI * 0.25;
                double[] u = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                        0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363 };
                double[] w = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                       0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
                for (int i = 0; i < 8; i++)
                {
                    double tau = half * (u[i] + 1.0);
                    double wt = half * w[i];
                    double s2v = Math.Pow(Math.Sin(tau), 2 * v);
                    var f = Complex.Sinh(Z * Math.Cos(tau));
                    sum += wt * s2v * f;
                }

                var pref = 2.0 * invNorm * Complex.Exp(v * Complex.Log(Z * 0.5));
                var val = pref * sum;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
            else
            {
                var lead = Complex.Exp((v - 1) * Complex.Log(Z * 0.5)) * invNorm;
                var iv = Special.I(x, v);
                var val = new Complex(iv.Real, iv.Imag) + lead;
                return new Complex32((float)val.Real, (float)val.Imaginary);
            }
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
            float r = 2.2360679774997896964f;
            float phi = (1.0f + r) / 2.0f;
            float psi = (1.0f - r) / 2.0f;
            float num = Maths.Pow(phi, n) - Maths.Pow(psi, n);
            return (int)(num / r);
        }
        /// <summary>
        /// Returns the value of the Luca number.
        /// </summary>
        /// <param name="n">Integer number</param>
        /// <returns>Integer number</returns>
        public static int Lucas(int n)
        {
            float r = 2.2360679774997896964f;
            float phi = (1.0f + r) / 2.0f;
            float psi = (1.0f - r) / 2.0f;
            float num = Maths.Pow(phi, n) + Maths.Pow(psi, n);
            return (int)(num);
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
            return Special.DiGamma(n) + 1.0f / n + Maths.Gamma;
        }
        /// <summary>
        /// Returns the harmonic number.
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="m">Value</param>
        /// <returns>Value</returns>
        public static float Harm(int n, float m)
        {
            float sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += Maths.Pow(i + 1, -m);
            }

            return sum;
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
            // properties:
            double p = 1, s = 0;
            double v = x - 0.5f;
            double u = Math.Pow(v, n);

            // series:
            for (int k = 0; k <= n; k++)
            {
                s += Special.Binomial(n, k) * Special.Euler(k) / p * u; p *= 2; u /= v;
            }

            // result:
            return s;
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
            // properties:
            double p = 1, s = 0;

            // series:
            for (int k = 0; k <= n; k++)
            {
                s += Special.Binomial(n, k) * Special.Bernoulli(n - k) * p; p *= x;
            }

            // result:
            return s;
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
            // Minkowski function:
            long p = x, q = 1, r = p + 1, s = 1, m, n;
            float d = 1, y = p;

            if (x < p || (p < 0) ^ (r <= 0))
            {
                return x;
            }

            // calculating:
            for (; ; )
            {
                d /= 2; if (y + d == y) break;
                m = p + r; if ((m < 0) ^ (p < 0)) break;
                n = q + s; if (n < 0) break;

                if (x < (float)m / n)
                {
                    r = m; s = n;
                }
                else
                {
                    y += d; p = m; q = n;
                }
            }
            return y + d;
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
