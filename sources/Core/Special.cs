using System;

namespace UMapx.Core
{
    /// <summary>
    /// Uses to implement special mathematical functions.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Special_functions
    /// </remarks>
    /// </summary>
    public static class Special
    {
        #region Private data
        /// <summary>
        /// 
        /// </summary>
        private const double LogMax = 7.09782712893383996732E2;
        /// <summary>
        /// 
        /// </summary>
        private const double LogMin = -7.451332191019412076235E2;
        /// <summary>
        /// 
        /// </summary>
        private const double sqrtPI = 1.7724538509055160272981674833411;
        #endregion

        #region Integral functions
        /// <summary>
        /// Returns the value of the integral cosine.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Ci(double x)
        {
            // special cases:
            if (x == 0)
                return double.NegativeInfinity;
            if (x < 0)
                return double.NaN;
            if (x > 35)
                return 0.0;

            // properties:
            double s = 0;
            double f = 1.0;
            double z = x * x;
            double t, m = 1, eps = 1e-16;
            int k, i, iterations = 120;
            int p = 1;

            // Taylor series:
            for (i = 1; i < iterations; i++)
            {
                // factorial:
                k = 2 * i;
                f *= k * (k - 1);

                // sign and value:
                p *= (-1);
                m *= z;
                t = p * m / (f * k);

                // stop point:
                if (Math.Abs(t) < eps)
                { break; }
                else { s += t; }
            }

            // construction:
            return Maths.Gamma + Math.Log(x) + s;
        }
        /// <summary>
        /// Returns the value of the integral sine.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Si(double x)
        {
            // special cases:
            if (x > 35)
                return 1.6;
            if (x < -35)
                return -1.6;

            // properties:
            double s = x;
            double f = 1.0;
            double z = x * x;
            double t, m = x, eps = 1e-16;
            int k, i, j = 1, iterations = 120;
            int p = 1;

            // Taylor series:
            for (i = 1; i < iterations; i++)
            {
                // factorial:
                k = 2 * i + 1; j += 2;
                f *= j * (k - 1);

                // sign and value:
                p *= (-1);
                m *= z;
                t = p / f / k * m;

                // stop point:
                if (Math.Abs(t) < eps)
                { break; }
                else { s += t; }
            }

            // result:
            return s;
        }
        /// <summary>
        /// Returns the value of an integral exponential function.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Ei(double x)
        {
            // Properties:
            double s = 0.0;
            double f = 1.0;
            double m = 1.0;
            double t, eps = 1e-8;
            int i, iterations = 120;

            // Taylor series:
            for (i = 1; i < iterations; i++)
            {
                // value:
                f *= i;
                m *= x;
                t = m / (f * i);

                // stop point:
                if (Math.Abs(t) < eps)
                { break; }
                else { s += t; }
            }

            // construction:
            double r = Maths.Gamma + s;

            // ranges:
            if (x < 0)
            {
                return r + Math.Log(-x);
            }
            return r + Math.Log(x);
        }
        /// <summary>
        /// Returns the value of the integral logarithm.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Li(double x)
        {
            // calculating Li(x) from Ei(x) 
            // integral function.

            if (x < 0)
            {
                return double.NaN;
            }
            return Ei(Maths.Log(x));
        }
        #endregion

        #region Fresnel integral functions
        /// <summary>
        /// Returns the value of the Fresnel integral C (x).
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Fresnelc(double x)
        {
            // special case:
            if (Math.Abs(x) > 6)
                return double.NaN;

            // properties:
            double s = x;
            double f = 1.0;
            double t, m = x, eps = 1e-16;
            double z = x * x * x * x;
            int k, i, j, iterations = 120;
            int p = 1;

            // Taylor series:
            for (i = 1; i < iterations; i++)
            {
                // factorial:
                j = 2 * i;
                k = 2 * j + 1;
                f *= j * (j - 1);

                // sign and value:
                p = -p;
                m *= z;
                t = p * m / (f * k);

                // stop point:
                if (Math.Abs(t) < eps)
                { break; }
                else { s += t; }
            }

            // construction:
            return s;
        }
        /// <summary>
        /// Returns the value of the Fresnel integral S(x).
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Fresnels(double x)
        {
            // special case:
            if (Math.Abs(x) > 6)
                return double.NaN;

            // properties:
            double s = x * x * x / 3;
            double f = 1.0;
            double t, m = 3.0 * s, eps = 1e-16;
            double z = x * x * x * x;
            int k, i, j = 1, iterations = 120;
            int p = 1;

            // Taylor series:
            for (i = 1; i < iterations; i++)
            {
                // factorial:
                k = 4 * i + 3;
                j = 2 * i;
                f *= j * (j + 1);

                // sign and value:
                p = -p;
                m *= z;
                t = p * m / f / k;

                // stop point:
                if (Math.Abs(t) < eps)
                { break; }
                else { s += t; }
            }

            // construction:
            return s;
        }
        #endregion

        #region Struve function
        /// <summary>
        /// Returns the value of the Struve function.
        /// </summary>
        /// <param name="x">Number</param>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double H(double x, int a)
        {
            // Struve function calculation method.
            // special cases:
            if (a < 0)
                return 0.0;

            // {1.0 / Г(3/2)} and {1.0 / Г(3/2 + a)}
            double x0 = 3.0 / 2.0, x1 = x0 + a;
            double g0 = 1.0 / Special.Gamma(x0);
            double g1 = 1.0 / Special.Gamma(x1);

            // construction:
            double b = x / 2.0;
            double s = 0, p = 1;
            double u = Math.Pow(b, a + 1);
            double t, eps = 1e-16;
            int i, iterations = 120;

            // Taylor series:
            for (i = 0; i < iterations; i++)
            {
                // value:
                t = p * g0 * g1 * u;

                // stop point:
                if (Math.Abs(t) < eps)
                {
                    break;
                }
                else
                {
                    // summary:
                    s += t;

                    // for next step:
                    g0 *= b / (x0 + i);
                    g1 *= b / (x1 + i);
                    p = -p;
                }
            }

            // result:
            return s;
        }
        /// <summary>
        /// Returns the value of the modified Struve function.
        /// </summary>
        /// <param name="x">Number</param>
        /// <param name="v">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double L(double x, int v)
        {
            // Modified Struve function calculation method.
            // special cases:
            if (v < 0)
                return 0.0;

            // {1.0 / Г(3/2)} and {1.0 / Г(3/2 + a)}
            double x0 = 3.0 / 2.0, x1 = x0 + v;
            double g0 = 1.0 / Special.Gamma(x0);
            double g1 = 1.0 / Special.Gamma(x1);

            // construction:
            double b = x / 2.0;
            double s = 0;
            double u = Math.Pow(b, v);
            double t, eps = 1e-16;
            int i, iterations = 120;

            // Taylor series:
            for (i = 0; i < iterations; i++)
            {
                // value:
                t = g0 * g1 * u;

                // stop point:
                if (Math.Abs(t) < eps)
                {
                    break;
                }
                else
                {
                    // summary:
                    s += t;

                    // for next step:
                    g0 *= b / (x0 + i);
                    g1 *= b / (x1 + i);
                }
            }

            // result:
            return b * u * s;
        }
        #endregion

        #region Beta function
        /// <summary>
        /// Returns the value of the beta function: B(a, b) = Г(a) * Г(b) / Г(ab).
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Beta(double a, double b)
        {
            return Special.Gamma(a) * Special.Gamma(b) / Special.Gamma(a + b);
        }
        /// <summary>
        /// Returns the value of the beta function: B(m, n) = (m - 1)! * (n - 1)! / (m + n - 1)!.
        /// </summary>
        /// <param name="m">Integer number</param>
        /// <param name="n">Integer number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Beta(int m, int n)
        {
            return Special.Factorial(m - 1) * Special.Factorial(n - 1) / Special.Factorial(m + n - 1);
        }
        /// <summary>
        /// Returns the value of an incomplete beta function: Bx(a, b).
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Number</param>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Beta(double a, double b, double x)
        {
            double aa, bb, t, xx, xc, w, y;
            bool flag;

            if (a <= 0.0)
                throw new ArgumentOutOfRangeException("Invalid argument value: a <= 0");
            if (b <= 0.0)
                throw new ArgumentOutOfRangeException("Invalid argument value: b <= 0");

            if ((x <= 0.0) || (x >= 1.0))
            {
                if (x == 0.0)
                {
                    return 0.0;
                }
                if (x == 1.0)
                {
                    return 1.0;
                }
                throw new ArgumentOutOfRangeException("Invalid argument value, beacause the value must belong to the interval [0, 1]");
            }

            flag = false;
            if ((b * x) <= 1.0 && x <= 0.95)
            {
                t = Special.Series(a, b, x);
                return t;
            }

            w = 1.0 - x;

            if (x > (a / (a + b)))
            {
                flag = true;
                aa = b;
                bb = a;
                xc = x;
                xx = w;
            }
            else
            {
                aa = a;
                bb = b;
                xc = w;
                xx = x;
            }

            if (flag && (bb * xx) <= 1.0 && xx <= 0.95)
            {
                t = Special.Series(aa, bb, xx);
                if (t <= double.Epsilon) t = 1.0 - double.Epsilon;
                else t = 1.0 - t;
                return t;
            }

            y = xx * (aa + bb - 2.0) - (aa - 1.0);
            if (y < 0.0)
                w = Special.Incbcf(aa, bb, xx);
            else
                w = Special.Incbd(aa, bb, xx) / xc;

            y = aa * System.Math.Log(xx);
            t = bb * System.Math.Log(xc);
            if ((aa + bb) < GammaMax && System.Math.Abs(y) < LogMax && System.Math.Abs(t) < LogMax)
            {
                t = System.Math.Pow(xc, bb);
                t *= System.Math.Pow(xx, aa);
                t /= aa;
                t *= w;
                t *= Special.Gamma(aa + bb) / (Special.Gamma(aa) * Special.Gamma(bb));
                if (flag)
                {
                    if (t <= double.Epsilon) t = 1.0 - double.Epsilon;
                    else t = 1.0 - t;
                }
                return t;
            }

            y += t + Special.GammaLog(aa + bb) - Special.GammaLog(aa) - Special.GammaLog(bb);
            y += System.Math.Log(w / aa);
            if (y < LogMin)
                t = 0.0;
            else
                t = System.Math.Exp(y);

            if (flag)
            {
                if (t <= double.Epsilon) t = 1.0 - double.Epsilon;
                else t = 1.0 - t;
            }
            return t;
        }
        /// <summary>
        /// Returns the value of a derivative beta function: B'(a, b).
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double BetaDerivative(double a, double b)
        {
            return Special.Beta(a, b) * (Special.DiGamma(a) - Special.DiGamma(a + b));
        }
        /// <summary>
        /// Returns the value of a regularized incomplete beta function: Ix(a, b).
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Number</param>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double BetaIncomplete(double a, double b, double x)
        {
            return Special.Beta(a, b, x) / Special.Beta(a, b);
        }
        #region Private beta approximations
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double Incbcf(double a, double b, double x)
        {
            double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
            double k1, k2, k3, k4, k5, k6, k7, k8;
            double r, t, ans, thresh;
            int n;
            double big = 4.503599627370496e15;
            double biginv = 2.22044604925031308085e-16;

            k1 = a;
            k2 = a + b;
            k3 = a;
            k4 = a + 1.0;
            k5 = 1.0;
            k6 = b - 1.0;
            k7 = k4;
            k8 = a + 2.0;

            pkm2 = 0.0;
            qkm2 = 1.0;
            pkm1 = 1.0;
            qkm1 = 1.0;
            ans = 1.0;
            r = 1.0;
            n = 0;
            thresh = 3.0 * double.Epsilon;

            do
            {
                xk = -(x * k1 * k2) / (k3 * k4);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                xk = (x * k5 * k6) / (k7 * k8);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                if (qk != 0) r = pk / qk;
                if (r != 0)
                {
                    t = System.Math.Abs((ans - r) / r);
                    ans = r;
                }
                else
                    t = 1.0;

                if (t < thresh) return ans;

                k1 += 1.0;
                k2 += 1.0;
                k3 += 2.0;
                k4 += 2.0;
                k5 += 1.0;
                k6 -= 1.0;
                k7 += 2.0;
                k8 += 2.0;

                if ((System.Math.Abs(qk) + System.Math.Abs(pk)) > big)
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }
                if ((System.Math.Abs(qk) < biginv) || (System.Math.Abs(pk) < biginv))
                {
                    pkm2 *= big;
                    pkm1 *= big;
                    qkm2 *= big;
                    qkm1 *= big;
                }
            } while (++n < 300);

            return ans;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double Incbd(double a, double b, double x)
        {
            double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
            double k1, k2, k3, k4, k5, k6, k7, k8;
            double r, t, ans, z, thresh;
            int n;
            double big = 4.503599627370496e15;
            double biginv = 2.22044604925031308085e-16;

            k1 = a;
            k2 = b - 1.0;
            k3 = a;
            k4 = a + 1.0;
            k5 = 1.0;
            k6 = a + b;
            k7 = a + 1.0;
            ;
            k8 = a + 2.0;

            pkm2 = 0.0;
            qkm2 = 1.0;
            pkm1 = 1.0;
            qkm1 = 1.0;
            z = x / (1.0 - x);
            ans = 1.0;
            r = 1.0;
            n = 0;
            thresh = 3.0 * double.Epsilon;
            do
            {
                xk = -(z * k1 * k2) / (k3 * k4);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                xk = (z * k5 * k6) / (k7 * k8);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                if (qk != 0) r = pk / qk;
                if (r != 0)
                {
                    t = System.Math.Abs((ans - r) / r);
                    ans = r;
                }
                else
                    t = 1.0;

                if (t < thresh) return ans;

                k1 += 1.0;
                k2 -= 1.0;
                k3 += 2.0;
                k4 += 2.0;
                k5 += 1.0;
                k6 += 1.0;
                k7 += 2.0;
                k8 += 2.0;

                if ((System.Math.Abs(qk) + System.Math.Abs(pk)) > big)
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }
                if ((System.Math.Abs(qk) < biginv) || (System.Math.Abs(pk) < biginv))
                {
                    pkm2 *= big;
                    pkm1 *= big;
                    qkm2 *= big;
                    qkm1 *= big;
                }
            } while (++n < 300);

            return ans;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double Series(double a, double b, double x)
        {
            double s, t, u, v, n, t1, z, ai;

            ai = 1.0 / a;
            u = (1.0 - b) * x;
            v = u / (a + 1.0);
            t1 = v;
            t = u;
            n = 2.0;
            s = 0.0;
            z = double.Epsilon * ai;
            while (System.Math.Abs(v) > z)
            {
                u = (n - b) * x / n;
                t *= u;
                v = t / (a + n);
                s += v;
                n += 1.0;
            }
            s += t1;
            s += ai;

            u = a * System.Math.Log(x);
            if ((a + b) < GammaMax && System.Math.Abs(u) < LogMax)
            {
                t = Gamma(a + b) / (Gamma(a) * Gamma(b));
                s = s * t * System.Math.Pow(x, a);
            }
            else
            {
                t = GammaLog(a + b) - GammaLog(a) - GammaLog(b) + u + System.Math.Log(s);
                if (t < LogMin) s = 0.0;
                else s = System.Math.Exp(t);
            }
            return s;
        }
        #endregion
        #endregion

        #region Gamma function
        /// <summary>
        /// 
        /// </summary>
        private const double GammaMax = 171.624376956302725;
        /// <summary>
        /// 
        /// </summary>
        private static double[] Px = {
                         1.60119522476751861407E-4,
                         1.19135147006586384913E-3,
                         1.04213797561761569935E-2,
                         4.76367800457137231464E-2,
                         2.07448227648435975150E-1,
                         4.94214826801497100753E-1,
                         9.99999999999999996796E-1
                     };
        /// <summary>
        /// 
        /// </summary>
        private static double[] Qx = {
                         -2.31581873324120129819E-5,
                         5.39605580493303397842E-4,
                         -4.45641913851797240494E-3,
                         1.18139785222060435552E-2,
                         3.58236398605498653373E-2,
                         -2.34591795718243348568E-1,
                         7.14304917030273074085E-2,
                         1.00000000000000000320E0
                     };
        /// <summary>
        /// Returns the value of the Euler Gamma function: Г(z).
        /// </summary>
        /// <param name="z">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gamma(double z)
        {
            double p, g;
            double q = Math.Abs(z);

            if (q > 33.0)
            {
                if (z < 0.0)
                {
                    p = Math.Floor(q);
                    if (p == q) return double.NaN;
                    g = q - p;
                    if (g > 0.5)
                    {
                        p += 1.0;
                        g = q - p;
                    }
                    g = q * Math.Sin(Math.PI * g);
                    if (g == 0.0) return double.NaN;
                    g = Math.Abs(g);
                    g = Math.PI / (g * Special.Stirling(q));

                    return -g;
                }
                else
                {
                    return Special.Stirling(z);
                }
            }

            g = 1.0;
            while (z >= 3.0)
            {
                z -= 1.0;
                g *= z;
            }

            while (z < 0.0)
            {
                if (z == 0.0)
                {
                    return double.NaN;
                }
                else if (z > -1.0E-9)
                {
                    return (g / ((1.0 + 0.5772156649015329 * z) * z));
                }
                g /= z;
                z += 1.0;
            }

            while (z < 2.0)
            {
                if (z == 0.0)
                {
                    return double.NaN;
                }
                else if (z < 1.0E-9)
                {
                    return (g / ((1.0 + 0.5772156649015329 * z) * z));
                }
                g /= z;
                z += 1.0;
            }

            if ((z == 2.0) || (z == 3.0)) return g;

            z -= 2.0;
            p = Special.Polynomials(z, Px, 6);
            q = Special.Polynomials(z, Qx, 7);
            return g * p / q;

        }
        /// <summary>
        /// Returns the value of the natural logarithm of the Euler Gamma function: ln[Г(z)].
        /// </summary>
        /// <param name="z">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double GammaLog(double z)
        {
            if (z < 0)
            {
                if (z > -6)
                {
                    return Math.Log(Gamma(z));
                }
                return double.NaN;
            }
            else if (z > 0)
            {
                return Special.GammaLogLanczos(z);
            }
            return double.NaN;
        }
        /// <summary>
        /// Returns the value of the Digamma function: ψ(z).
        /// </summary>
        /// <param name="z">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double DiGamma(double z)
        {
            if (z == 0)
                return Double.NegativeInfinity;

            double s = 0;
            double w = 0;
            double y = 0;
            double yy = 0;
            double nz = 0;

            bool negative = false;

            if (z <= 0.0)
            {
                negative = true;
                double q = z;
                double p = (int)System.Math.Floor(q);

                if (p == q)
                    return double.NaN;

                nz = q - p;

                if (nz != 0.5)
                {
                    if (nz > 0.5)
                    {
                        p = p + 1.0;
                        nz = q - p;
                    }
                    nz = Math.PI / Math.Tan(System.Math.PI * nz);
                }
                else
                {
                    nz = 0.0;
                }

                z = 1.0 - z;
            }

            if (z <= 10.0 & z == Math.Floor(z))
            {
                y = 0.0;
                int n = (int)Math.Floor(z);
                for (int i = 1; i <= n - 1; i++)
                {
                    w = i;
                    y = y + 1.0 / w;
                }
                y = y - 0.57721566490153286061;
            }
            else
            {
                s = z;
                w = 0.0;

                while (s < 10.0)
                {
                    w = w + 1.0 / s;
                    s = s + 1.0;
                }

                if (s < 1.0E17)
                {
                    yy = 1.0 / (s * s);

                    double polv = 8.33333333333333333333E-2;
                    polv = polv * yy - 2.10927960927960927961E-2;
                    polv = polv * yy + 7.57575757575757575758E-3;
                    polv = polv * yy - 4.16666666666666666667E-3;
                    polv = polv * yy + 3.96825396825396825397E-3;
                    polv = polv * yy - 8.33333333333333333333E-3;
                    polv = polv * yy + 8.33333333333333333333E-2;
                    y = yy * polv;
                }
                else
                {
                    y = 0.0;
                }
                y = Math.Log(s) - 0.5 / s - y - w;
            }

            if (negative == true)
            {
                y = y - nz;
            }

            return y;
        }
        /// <summary>
        /// Returns the value of the Trigamma function: ψ1(z).
        /// </summary>
        /// <param name="z">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double TriGamma(double z)
        {
            double a = 0.0001;
            double b = 5.0;
            double b2 = 0.1666666667;
            double b4 = -0.03333333333;
            double b6 = 0.02380952381;
            double b8 = -0.03333333333;
            double value;
            double y;
            double yy;

            // Check the input.
            if (z <= 0.0)
            {
                return double.NaN;
            }

            yy = z;

            // Use small value approximation if X <= A.
            if (z <= a)
            {
                value = 1.0 / z / z;
                return value;
            }

            // Increase argument to ( X + I ) >= B.
            value = 0.0;

            while (yy < b)
            {
                value = value + 1.0 / yy / yy;
                yy = yy + 1.0;
            }

            // Apply asymptotic formula if argument is B or greater.
            y = 1.0 / yy / yy;

            value = value + 0.5 *
                y + (1.0
              + y * (b2
              + y * (b4
              + y * (b6
              + y * b8)))) / yy;

            return value;
        }
        /// <summary>
        /// Returns the value of the degree of the Euler Gamma function: Г(z)^p.
        /// </summary>
        /// <param name="z">Number</param>
        /// <param name="p">Power</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gamma(double z, uint p)
        {
            if (p == 1)
            {
                return Special.Gamma(z);
            }

            double prod = Math.Pow(Math.PI, (1 / 4.0) * p * (p - 1));
            int i;

            for (i = 0; i < p; i++)
            {
                prod *= Special.Gamma(z - 0.5 * i);
            }

            return prod;
        }
        /// <summary>
        /// Returns the value of the incomplete upper Gamma function: Q(s, x) = Γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Number</param>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double GammaQ(double s, double x)
        {
            const double big = 4.503599627370496e15;
            const double biginv = 2.22044604925031308085e-16;
            double ans, ax, c, yc, r, t, y, z;
            double pk, pkm1, pkm2, qk, qkm1, qkm2;

            if (x <= 0 || s <= 0)
                return 1.0;

            if (x < 1.0 || x < s)
                return 1.0 - GammaP(s, x);

            if (Double.IsPositiveInfinity(x))
                return 0;

            ax = s * Math.Log(x) - x - GammaLog(s);

            if (ax < -LogMax)
                return 0.0;

            ax = Math.Exp(ax);

            // continued fraction
            y = 1.0 - s;
            z = x + y + 1.0;
            c = 0.0;
            pkm2 = 1.0;
            qkm2 = x;
            pkm1 = x + 1.0;
            qkm1 = z * x;
            ans = pkm1 / qkm1;

            do
            {
                c += 1.0;
                y += 1.0;
                z += 2.0;
                yc = y * c;
                pk = pkm1 * z - pkm2 * yc;
                qk = qkm1 * z - qkm2 * yc;
                if (qk != 0)
                {
                    r = pk / qk;
                    t = Math.Abs((ans - r) / r);
                    ans = r;
                }
                else
                    t = 1.0;

                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if (Math.Abs(pk) > big)
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }
            } while (t > double.Epsilon);

            return ans * ax;
        }
        /// <summary>
        /// Returns the value of an incomplete lower Gamma function: P(s, x) = γ(s, x) / Γ(s).
        /// </summary>
        /// <param name="s">Number</param>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double GammaP(double s, double x)
        {
            if (s <= 0)
                return 1.0;

            if (x <= 0)
                return 0.0;

            if (x > 1.0 && x > s)
                return 1.0 - GammaQ(s, x);

            double ax = s * Math.Log(x) - x - GammaLog(s);

            if (ax < -LogMax)
                return 0.0;

            ax = Math.Exp(ax);

            double r = s;
            double c = 1.0;
            double ans = 1.0;

            do
            {
                r += 1.0;
                c *= x / r;
                ans += c;
            } while (c / ans > double.Epsilon);

            return ans * ax / s;
        }
        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x).
        /// </summary>
        /// <param name="s">Number</param>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double GammaIncomplete(double s, double x)
        {
            double ans, ax, c, r;

            if (x <= 0 || s <= 0) return 0.0;

            if (x > 1.0 && x > s) return 1.0 - GammaIncomplete(s, x, true);

            /* Compute  x**a * exp(-x) / gamma(a)  */
            ax = s * Math.Log(x) - x - GammaLog(s);
            if (ax < -LogMax) return (0.0);

            ax = Math.Exp(ax);

            /* power series */
            r = s;
            c = 1.0;
            ans = 1.0;

            do
            {
                r += 1.0;
                c *= x / r;
                ans += c;
            } while (c / ans > double.Epsilon);

            return (ans * ax / s);

        }
        /// <summary>
        /// Returns the value of an incomplete Gamma function: γ(s, x).
        /// </summary>
        /// <param name="s">Number</param>
        /// <param name="x">Number</param>
        /// <param name="complemented">Additional function or not</param>
        /// <returns>Double precision floating point number</returns>
        public static double GammaIncomplete(double s, double x, bool complemented)
        {
            // not complemented
            if (!complemented)
                return GammaIncomplete(s, x);

            double big = 4.503599627370496e15;
            double biginv = 2.22044604925031308085e-16;
            double ans, ax, c, yc, r, t, y, z;
            double pk, pkm1, pkm2, qk, qkm1, qkm2;

            if (x <= 0 || s <= 0) return 1.0;

            if (x < 1.0 || x < s) return 1.0 - GammaIncomplete(s, x);

            ax = s * Math.Log(x) - x - GammaLog(s);
            if (ax < -LogMax) return 0.0;

            ax = Math.Exp(ax);

            /* continued fraction */
            y = 1.0 - s;
            z = x + y + 1.0;
            c = 0.0;
            pkm2 = 1.0;
            qkm2 = x;
            pkm1 = x + 1.0;
            qkm1 = z * x;
            ans = pkm1 / qkm1;

            do
            {
                c += 1.0;
                y += 1.0;
                z += 2.0;
                yc = y * c;
                pk = pkm1 * z - pkm2 * yc;
                qk = qkm1 * z - qkm2 * yc;
                if (qk != 0)
                {
                    r = pk / qk;
                    t = Math.Abs((ans - r) / r);
                    ans = r;
                }
                else
                    t = 1.0;

                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if (Math.Abs(pk) > big)
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }
            } while (t > double.Epsilon);

            return ans * ax;
        }

        #region Gamma approximations
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="coef"></param>
        /// <param name="N"></param>
        /// <returns></returns>
        private static double Polynomials(double x, double[] coef, int N)
        {
            double sum = coef[0];

            for (int i = 1; i <= N; i++)
            {
                sum = sum * x + coef[i];
            }

            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double Stirling(double x)
        {
            double[] STIR =
            {
                7.87311395793093628397E-4,
               -2.29549961613378126380E-4,
               -2.68132617805781232825E-3,
                3.47222221605458667310E-3,
                8.33333333333482257126E-2,
            };
            double MAXSTIR = 143.01608;

            double w = 1.0 / x;
            double y = Math.Exp(x);

            w = 1.0 + w * Polynomials(w, STIR, 4);

            if (x > MAXSTIR)
            {
                /* Avoid overflow in Math.Pow() */
                double v = Math.Pow(x, 0.5 * x - 0.25);
                y = v * (v / y);
            }
            else
            {
                y = Math.Pow(x, x - 0.5) / y;
            }
            y = 2.50662827463100050242E0 * y * w;
            return y;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="z"></param>
        /// <returns></returns>
        private static double GammaLogLanczos(double z)
        {
            double[] coef = new double[6] { 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179E-2, -0.5395239384953E-5 };
            double LogSqrtTwoPi = 0.91893853320467274178;
            double denom = z + 1;
            double y = z + 5.5;
            double series = 1.000000000190015;
            int i;

            for (i = 0; i < 6; ++i)
            {
                series += coef[i] / denom;
                denom += 1.0;
            }
            return (LogSqrtTwoPi + (z + 0.5) * Math.Log(y) -
            y + Math.Log(series / z));
        }
        #endregion
        #endregion

        #region Laplace function
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erf(double x)
        {
            return Erf(x, false);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="inverse">Reverse function or not</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erf(double x, bool inverse)
        {
            // exception
            if (x == 0.0)
                return 0.0;

            // params
            int s = inverse ? 1 : -1;
            double y = s * x * x;
            double z = 2 / Math.Sqrt(Math.PI);
            double a = 1, b = x, c;
            double eps = 1e-16;
            int i, iterations = 120, k = 3;

            // Taylor series:
            for (i = 1; i < iterations; i++, k += 2)
            {
                // series:
                a *= y / i;
                c = a * x / k;

                // stop point:
                if (Math.Abs(c) < eps)
                {
                    break;
                }
                else
                {
                    b += c;
                }
            }

            // erf(x)
            double erf = z * b;
            if (!inverse)
            {
                if (Math.Abs(erf) > 1.0)
                    return Math.Sign(erf);
            }

            return erf;
        }
        /// <summary>
        /// Returns the value of the imaginary error function.
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erfi(double x)
        {
            // special cases:
            if (x > 26.658987628)
                return double.PositiveInfinity;
            if (x < -26.658987628)
                return double.NegativeInfinity;

            // properties:
            double s = x;
            double m = x, t;
            double z = x * x;
            int k, iterations = 930;
            double eps = 1e-16;
            double p = 2.0 / Special.sqrtPI;

            // Taylor series:
            for (int i = 1; i < iterations; i++)
            {
                // value:
                k = 2 * i + 1;
                m *= z / i;
                t = m / k;

                // stop point:
                if (Math.Abs(t) < eps)
                { break; }
                else { s += t; }
            }

            // construction:
            return p * s;
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="a">The lower boundary of the normalization</param>
        /// <param name="b">The upper limit of the normalization</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erf(double x, double a, double b)
        {
            return Erf((x - a) / b);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="range">A pair of fractional numbers</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erf(double x, RangeDouble range)
        {
            return Erf(x, range.Min, range.Max);
        }

        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erfc(double x)
        {
            return 1.0 - Erf(x);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="a">The lower boundary of the normalization</param>
        /// <param name="b">The upper limit of the normalization</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erfc(double x, double a, double b)
        {
            return 1.0 - Erf(x, a, b);
        }
        /// <summary>
        /// Returns the value of the Laplace integral (an additional error function).
        /// </summary>
        /// <param name="x">The value of the upper limit of the integral</param>
        /// <param name="range">A pair of fractional numbers</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erfc(double x, RangeDouble range)
        {
            return 1.0 - Erf(x, range.Min, range.Max);
        }
        #endregion

        #region Q-function
        /// <summary>
        /// Returns the value of a Q function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="inverse">Inverse function or not</param>
        /// <returns>Double precision floating point number</returns>
        public static double Q(double x, bool inverse = false)
        {
            if (inverse)
            {
                return Math.Sqrt(2) * Special.Erf(1 - 2 * x, true);
            }
            return 0.5 * Special.Erfc(x / Math.Sqrt(2));
        }
        #endregion

        #region Bessel functions
        /// <summary>
        /// Returns the value of the Bessel function of the first kind at a = 0.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        private static double J0(double x)
        {
            double ax;

            if ((ax = System.Math.Abs(x)) < 8.0)
            {
                double y = x * x;
                double ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
                    + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
                double ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                    + y * (59272.64853 + y * (267.8532712 + y * 1.0))));

                return ans1 / ans2;
            }
            else
            {
                double z = 8.0 / ax;
                double y = z * z;
                double xx = ax - 0.785398164;
                double ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                    + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                double ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                    + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                    - y * 0.934935152e-7)));

                return System.Math.Sqrt(0.636619772 / ax) *
                    (System.Math.Cos(xx) * ans1 - z * System.Math.Sin(xx) * ans2);
            }
        }
        /// <summary>
        /// Returns the value of the Bessel function of the first kind at a = 1.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        private static double J1(double x)
        {
            double ax;
            double y;
            double ans1, ans2;

            if ((ax = System.Math.Abs(x)) < 8.0)
            {
                y = x * x;
                ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
                    + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
                ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
                    + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
                return ans1 / ans2;
            }
            else
            {
                double z = 8.0 / ax;
                double xx = ax - 2.356194491;
                y = z * z;

                ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                    + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                ans2 = 0.04687499995 + y * (-0.2002690873e-3
                    + y * (0.8449199096e-5 + y * (-0.88228987e-6
                    + y * 0.105787412e-6)));
                double ans = System.Math.Sqrt(0.636619772 / ax) *
                    (System.Math.Cos(xx) * ans1 - z * System.Math.Sin(xx) * ans2);
                if (x < 0.0) ans = -ans;
                return ans;
            }
        }
        /// <summary>
        /// Returns the value of a Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Number</param>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double J(double x, int a)
        {
            if (a < 0 || x == 0) return 0;
            if (a == 0) return J0(x);
            if (a == 1) return J1(x);

            int j, m;
            double ax, bj, bjm, bjp, sum, tox, ans;
            bool jsum;
            double ACC = 40.0;
            double BIGNO = 1.0e+10;
            double BIGNI = 1.0e-10;

            ax = System.Math.Abs(x);
            if (ax == 0.0) return 0.0;
            else if (ax > (double)a)
            {
                tox = 2.0 / ax;
                bjm = J0(ax);
                bj = J1(ax);
                for (j = 1; j < a; j++)
                {
                    bjp = j * tox * bj - bjm;
                    bjm = bj;
                    bj = bjp;
                }
                ans = bj;
            }
            else
            {
                tox = 2.0 / ax;
                m = 2 * ((a + (int)System.Math.Sqrt(ACC * a)) / 2);
                jsum = false;
                bjp = ans = sum = 0.0;
                bj = 1.0;
                for (j = m; j > 0; j--)
                {
                    bjm = j * tox * bj - bjp;
                    bjp = bj;
                    bj = bjm;
                    if (System.Math.Abs(bj) > BIGNO)
                    {
                        bj *= BIGNI;
                        bjp *= BIGNI;
                        ans *= BIGNI;
                        sum *= BIGNI;
                    }
                    if (jsum) sum += bj;
                    jsum = !jsum;
                    if (j == a) ans = bjp;
                }
                sum = 2.0 * sum - bj;
                ans /= sum;
            }

            return x < 0.0 && a % 2 == 1 ? -ans : ans;
        }
        /// <summary>
        /// Returns the value of the Bessel function of the second kind at a = 0.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        private static double Y0(double x)
        {
            if (x < 8.0)
            {
                double y = x * x;

                double ans1 = -2957821389.0 + y * (7062834065.0 + y * (-512359803.6
                    + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));
                double ans2 = 40076544269.0 + y * (745249964.8 + y * (7189466.438
                    + y * (47447.26470 + y * (226.1030244 + y * 1.0))));

                return (ans1 / ans2) + 0.636619772 * J0(x) * System.Math.Log(x);
            }
            else
            {
                double z = 8.0 / x;
                double y = z * z;
                double xx = x - 0.785398164;

                double ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                    + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                double ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                    + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                    + y * (-0.934945152e-7))));
                return System.Math.Sqrt(0.636619772 / x) *
                    (System.Math.Sin(xx) * ans1 + z * System.Math.Cos(xx) * ans2);
            }
        }
        /// <summary>
        /// Returns the value of the Bessel function of the second kind at a = 1.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        private static double Y1(double x)
        {
            if (x < 8.0)
            {
                double y = x * x;
                double ans1 = x * (-0.4900604943e13 + y * (0.1275274390e13
                    + y * (-0.5153438139e11 + y * (0.7349264551e9
                    + y * (-0.4237922726e7 + y * 0.8511937935e4)))));
                double ans2 = 0.2499580570e14 + y * (0.4244419664e12
                    + y * (0.3733650367e10 + y * (0.2245904002e8
                    + y * (0.1020426050e6 + y * (0.3549632885e3 + y)))));
                return (ans1 / ans2) + 0.636619772 * (J1(x) * System.Math.Log(x) - 1.0 / x);
            }
            else
            {
                double z = 8.0 / x;
                double y = z * z;
                double xx = x - 2.356194491;
                double ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
                    + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
                double ans2 = 0.04687499995 + y * (-0.2002690873e-3
                    + y * (0.8449199096e-5 + y * (-0.88228987e-6
                    + y * 0.105787412e-6)));
                return System.Math.Sqrt(0.636619772 / x) *
                    (System.Math.Sin(xx) * ans1 + z * System.Math.Cos(xx) * ans2);
            }
        }
        /// <summary>
        /// Returns the value of a Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Number</param>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Y(double x, int a)
        {
            if (a < 0 || x == 0) return 0;
            if (a == 0) return Y0(x);
            if (a == 1) return Y1(x);

            double by, bym, byp, tox;
            tox = 2.0 / x;
            by = Y1(x);
            bym = Y0(x);
            for (int j = 1; j < a; j++)
            {
                byp = j * tox * by - bym;
                bym = by;
                by = byp;
            }
            return by;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the first kind at a = 0.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        private static double I0(double x)
        {
            double ans;
            double ax = Math.Abs(x);

            if (ax < 3.75)
            {
                double y = x / 3.75;
                y = y * y;
                ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                   + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
            }
            else
            {
                double y = 3.75 / ax;
                ans = (Math.Exp(ax) / Math.Sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
                   + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
                   + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
                   + y * 0.392377e-2))))))));
            }

            return ans;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the first kind at a = 1.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        private static double I1(double x)
        {
            double ans;

            double ax = Math.Abs(x);

            if (ax < 3.75)
            {
                double y = x / 3.75;
                y = y * y;
                ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
                   + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
            }
            else
            {
                double y = 3.75 / ax;
                ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
                ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
                ans *= Math.Exp(ax) / Math.Sqrt(ax);
            }
            return x < 0.0 ? -ans : ans;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="x">Number</param>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double I(double x, int a)
        {
            if (a < 0 || x == 0) return 0;
            if (a == 0) return I0(x);
            if (a == 1) return I1(x);

            double ACC = 40.0;
            double BIGNO = 1.0e+10;
            double BIGNI = 1.0e-10;
            double tox = 2.0 / Math.Abs(x);
            double bip = 0, ans = 0.0;
            double bi = 1.0;

            for (int j = 2 * (a + (int)Math.Sqrt(ACC * a)); j > 0; j--)
            {
                double bim = bip + j * tox * bi;
                bip = bi;
                bi = bim;

                if (Math.Abs(bi) > BIGNO)
                {
                    ans *= BIGNI;
                    bi *= BIGNI;
                    bip *= BIGNI;
                }

                if (j == a)
                    ans = bip;
            }

            ans *= I0(x) / bi;
            return x < 0.0 && a % 2 == 1 ? -ans : ans;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the second kind at a = 0.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        private static double K0(double x)
        {
            double y, ans;

            if (x <= 2.0)
            {
                y = x * x / 4.0;
                ans = (-Math.Log(x / 2.0) * I0(x)) + (-0.57721566 + y * (0.42278420
                   + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2
                   + y * (0.10750e-3 + y * 0.74e-5))))));
            }
            else
            {
                y = 2.0 / x;
                ans = (Math.Exp(-x) / Math.Sqrt(x)) * (1.25331414 + y * (-0.7832358e-1
                   + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2
                   + y * (-0.251540e-2 + y * 0.53208e-3))))));
            }
            return ans;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the second kind at a = 1.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        private static double K1(double x)
        {
            double y, ans;

            if (x <= 2.0)
            {
                y = x * x / 4.0;
                ans = (Math.Log(x / 2.0) * I1(x)) + (1.0 / x) * (1.0 + y * (0.15443144
                   + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1
                   + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
            }
            else
            {
                y = 2.0 / x;
                ans = (Math.Exp(-x) / Math.Sqrt(x)) * (1.25331414 + y * (0.23498619
                   + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2
                   + y * (0.325614e-2 + y * (-0.68245e-3)))))));
            }
            return ans;
        }
        /// <summary>
        /// Returns the value of the modified Bessel function of the second kind.
        /// </summary>
        /// <param name="x">Number</param>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double K(double x, int a)
        {
            if (a < 0 || x == 0) return 0;
            if (a == 0) return (K0(x));
            if (a == 1) return (K1(x));

            int j;
            double bk, bkm, bkp, tox;
            tox = 2.0 / x;
            bkm = K0(x);
            bk = K1(x);
            for (j = 1; j < a; j++)
            {
                bkp = bkm + j * tox * bk;
                bkm = bk;
                bk = bkp;
            }
            return bk;
        }
        #endregion

        #region Owen's T-function components
        /// <summary>
        /// Returns the value of the Owen T function.
        /// </summary>
        /// <param name="h">First argument</param>
        /// <param name="a">Second argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Owen(double h, double a)
        {
            double absa;
            double absh;
            double ah;
            double cut = 0.67;
            double normah;
            double normh;
            double value;

            absh = Math.Abs(h);
            absa = Math.Abs(a);
            ah = absa * absh;

            if (absa <= 1.0)
            {
                value = Special.owenhaah(absh, absa, ah);
            }
            else if (absh <= cut)
            {
                value = 0.25 - (-0.5 + Special.Q(-absh)) * (-0.5 + Special.Q(-ah))
                  - Special.owenhaah(ah, 1.0 / absa, absh);
            }
            else
            {
                normh = Special.Q(absh);
                normah = Special.Q(ah);
                value = 0.5 * (normh + normah) - normh * normah
                - Special.owenhaah(ah, 1.0 / absa, absh);
            }

            if (a < 0.0)
            {
                value = -value;
            }

            return value;
        }
        #region Private data
        /// <summary>
        /// Returns the value of the Owen T function.
        /// </summary>
        /// <param name="h">First argument</param>
        /// <param name="a">Second argument</param>
        /// <param name="ah">h * a</param>
        /// <returns>Double precision floating point number</returns>
        private static double owenhaah(double h, double a, double ah)
        {
            double ai;
            double aj;
            double AS;
            double dhs;
            double dj;
            double gj;

            double hs;
            int i;
            int iaint;
            int icode;
            int ihint;
            int ii;
            int j;
            int jj;
            int m;
            int maxii;
            double normh;

            double r;
            const double rrtpi = 0.39894228040143267794;
            const double rtwopi = 0.15915494309189533577;

            double y;
            double yi;
            double z;
            double zi;

            double value = 0;
            double vi;


            /*
              Determine appropriate method from t1...t6
            */
            ihint = 15;

            for (i = 1; i <= 14; i++)
            {
                if (h <= hrange[i - 1])
                {
                    ihint = i;
                    break;
                }
            }

            iaint = 8;

            for (i = 1; i <= 7; i++)
            {
                if (a <= arange[i - 1])
                {
                    iaint = i;
                    break;
                }
            }

            icode = select[ihint - 1 + (iaint - 1) * 15];
            m = ord[icode - 1];

            /*
              t1(h, a, m) ; m = 2, 3, 4, 5, 7, 10, 12 or 18
              jj = 2j - 1 ; gj = exp(-h*h/2) * (-h*h/2)**j / j
              aj = a**(2j-1) / (2*pi)
            */

            if (meth[icode - 1] == 1)
            {
                hs = -0.5 * h * h;
                dhs = Math.Exp(hs);
                AS = a * a;
                j = 1;
                jj = 1;
                aj = rtwopi * a;
                value = rtwopi * Math.Atan(a);
                dj = dhs - 1.0;
                gj = hs * dhs;

                for (; ; )
                {
                    value = value + dj * aj / (double)(jj);

                    if (m <= j)
                    {
                        return value;
                    }
                    j = j + 1;
                    jj = jj + 2;
                    aj = aj * AS;
                    dj = gj - dj;
                    gj = gj * hs / (double)(j);
                }
            }

            /*
              t2(h, a, m) ; m = 10, 20 or 30
              z = (-1)^(i-1) * zi ; ii = 2i - 1
              vi = (-1)^(i-1) * a^(2i-1) * exp[-(a*h)^2/2] / sqrt(2*pi)
            */
            else if (meth[icode - 1] == 2)
            {
                maxii = m + m + 1;
                ii = 1;
                value = 0.0;
                hs = h * h;
                AS = -a * a;
                vi = rrtpi * a * Math.Exp(-0.5 * ah * ah);
                z = 0.5 * (-0.5 + Special.Q(-ah)) / h;
                y = 1.0 / hs;

                for (; ; )
                {
                    value = value + z;

                    if (maxii <= ii)
                    {
                        value = value * rrtpi * Math.Exp(-0.5 * hs);
                        return value;
                    }
                    z = y * (vi - (double)(ii) * z);
                    vi = AS * vi;
                    ii = ii + 2;
                }
            }
            /*
              t3(h, a, m) ; m = 20
              ii = 2i - 1
              vi = a**(2i-1) * exp[-(a*h)**2/2] / sqrt(2*pi)
            */
            else if (meth[icode - 1] == 3)
            {
                i = 1;
                ii = 1;
                value = 0.0;
                hs = h * h;
                AS = a * a;
                vi = rrtpi * a * Math.Exp(-0.5 * ah * ah);
                zi = 0.5 * (-0.5 + Special.Q(-ah)) / h;
                y = 1.0 / hs;

                for (; ; )
                {
                    value = value + zi * coefT[i - 1];

                    if (m < i)
                    {
                        value = value * rrtpi * Math.Exp(-0.5 * hs);
                        return value;
                    }
                    zi = y * ((double)(ii) * zi - vi);
                    vi = AS * vi;
                    i = i + 1;
                    ii = ii + 2;
                }
            }
            /*
              t4(h, a, m) ; m = 4, 7, 8 or 20;  ii = 2i + 1
              ai = a * exp[-h*h*(1+a*a)/2] * (-a*a)**i / (2*pi)
            */
            else if (meth[icode - 1] == 4)
            {
                maxii = m + m + 1;
                ii = 1;
                hs = h * h;
                AS = -a * a;
                value = 0.0;
                ai = rtwopi * a * Math.Exp(-0.5 * hs * (1.0 - AS));
                yi = 1.0;

                for (; ; )
                {
                    value = value + ai * yi;

                    if (maxii <= ii)
                        return value;

                    ii = ii + 2;
                    yi = (1.0 - hs * yi) / (double)(ii);
                    ai = ai * AS;
                }
            }
            /*
              t5(h, a, m) ; m = 13
              2m - point gaussian quadrature
            */
            else if (meth[icode - 1] == 5)
            {
                value = 0.0;
                AS = a * a;
                hs = -0.5 * h * h;
                for (i = 1; i <= m; i++)
                {
                    r = 1.0 + AS * pts[i - 1];
                    value = value + wts[i - 1] * Math.Exp(hs * r) / r;
                }
                value = a * value;
            }
            /*
              t6(h, a);  approximation for a near 1, (a<=1)
            */
            else if (meth[icode - 1] == 6)
            {
                normh = Special.Q(h);
                value = 0.5 * normh * (1.0 - normh);
                y = 1.0 - a;
                r = Math.Atan(y / (1.0 + a));

                if (r != 0.0)
                    value = value - rtwopi * r * Math.Exp(-0.5 * y * h * h / r);
            }

            return value;
        }
        /// <summary>
        /// 
        /// </summary>
        private static readonly double[] arange =
        {
            0.025, 0.09, 0.15, 0.36, 0.5, 0.9, 0.99999
        };
        /// <summary>
        /// 
        /// </summary>
        private static readonly double[] coefT =
        {
                                          0.99999999999999987510,
            -0.99999999999988796462,      0.99999999998290743652,
            -0.99999999896282500134,      0.99999996660459362918,
            -0.99999933986272476760,      0.99999125611136965852,
            -0.99991777624463387686,      0.99942835555870132569,
            -0.99697311720723000295,      0.98751448037275303682,
            -0.95915857980572882813,      0.89246305511006708555,
            -0.76893425990463999675,      0.58893528468484693250,
            -0.38380345160440256652,      0.20317601701045299653,
            -0.82813631607004984866E-01,  0.24167984735759576523E-01,
            -0.44676566663971825242E-02,  0.39141169402373836468E-03
        };
        /// <summary>
        /// 
        /// </summary>
        private static readonly double[] hrange =
        {
            0.02, 0.06, 0.09, 0.125, 0.26,
            0.4,  0.6,  1.6,  1.7,   2.33,
            2.4,  3.36, 3.4,  4.8
        };
        /// <summary>
        /// 
        /// </summary>
        private static readonly int[] meth =
        {
            1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6
        };
        /// <summary>
        /// 
        /// </summary>
        private static readonly int[] ord =
        {
            2, 3, 4, 5, 7,10,12,18,10,20,30,20, 4, 7, 8,20,13, 0
        };
        /// <summary>
        /// 
        /// </summary>
        private static readonly double[] pts =
        {
                                         0.35082039676451715489E-02,
            0.31279042338030753740E-01,  0.85266826283219451090E-01,
            0.16245071730812277011,      0.25851196049125434828,
            0.36807553840697533536,      0.48501092905604697475,
            0.60277514152618576821,      0.71477884217753226516,
            0.81475510988760098605,      0.89711029755948965867,
            0.95723808085944261843,      0.99178832974629703586
        };
        /// <summary>
        /// 
        /// </summary>
        private static readonly int[] select =
        {
            1, 1, 2,13,13,13,13,13,13,13,13,16,16,16, 9,
            1, 2, 2, 3, 3, 5, 5,14,14,15,15,16,16,16, 9,
            2, 2, 3, 3, 3, 5, 5,15,15,15,15,16,16,16,10,
            2, 2, 3, 5, 5, 5, 5, 7, 7,16,16,16,16,16,10,
            2, 3, 3, 5, 5, 6, 6, 8, 8,17,17,17,12,12,11,
            2, 3, 5, 5, 5, 6, 6, 8, 8,17,17,17,12,12,12,
            2, 3, 4, 4, 6, 6, 8, 8,17,17,17,17,17,12,12,
            2, 3, 4, 4, 6, 6,18,18,18,18,17,17,17,12,12
        };
        /// <summary>
        /// 
        /// </summary>
        private static readonly double[] wts =
        {
                                         0.18831438115323502887E-01,
            0.18567086243977649478E-01,  0.18042093461223385584E-01,
            0.17263829606398753364E-01,  0.16243219975989856730E-01,
            0.14994592034116704829E-01,  0.13535474469662088392E-01,
            0.11886351605820165233E-01,  0.10070377242777431897E-01,
            0.81130545742299586629E-02,  0.60419009528470238773E-02,
            0.38862217010742057883E-02,  0.16793031084546090448E-02
        };
        #endregion
        #endregion

        #region Guderman & Hartley functions
        /// <summary>
        /// Returns the value of the inverse Guderman function.
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Agd(double x)
        {
            return 0.5 * Math.Log((1.0 + Math.Sin(x)) / (1.0 - Math.Sin(x)));
        }
        /// <summary>
        /// Returns the value of the Guderman function.
        /// </summary>
        /// <param name="x">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gd(double x)
        {
            return 2 * Maths.Atg(Maths.Exp(x)) - 1.57079632679f;
        }
        /// <summary>
        /// Returns the value of the function Cas(x).
        /// </summary>
        /// <param name="theta">Theta</param>
        /// <returns>Double precision floating point number</returns>
        public static double Cas(double theta)
        {
            return 1.414f * Math.Sin(theta + 0.78539);
        }
        #endregion

        #region Sinc function
        /// <summary>
        /// Returns the value of the normalized cardinal sine function: f(x) = sin(πx) / (πx).
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sinc(double x)
        {
            return Special.Sinc(x, Math.PI);
        }
        /// <summary>
        /// Returns the value of the cardinal sine function with the parameter: f(x, a) = sin(ax) / (ax).
        /// </summary>
        /// <param name="x">Number</param>
        /// <param name="a">Parameter</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sinc(double x, double a)
        {
            if (x == 0)
            {
                return 1;
            }
            double ax = a * x;
            return Math.Sin(ax) / (ax);
        }
        #endregion

        #region Binomial function
        /// <summary>
        /// Returns the value of binomial coefficients: C(n, k) = n! / k! / (n-k)! для k > 0.
        /// </summary>
        /// <param name="n">Number</param>
        /// <param name="k">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Binomial(double n, double k)
        {
            if (k < 0)
            {
                return 0;
            }
            return Special.Factorial(n) / Special.Factorial(k) / Special.Factorial(n - k);
        }
        /// <summary>
        /// Returns the natural logarithm of binomial coefficients: log(C(n, k)) = log(n!) - log(k!) - log(n-k!).
        /// </summary>
        /// <param name="n">Number</param>
        /// <param name="k">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double LogBinomial(double n, double k)
        {
            return Special.LogFactorial(n) - Special.LogFactorial(k) - Special.LogFactorial(n - k);
        }
        #endregion

        #region Factorial function
        /// <summary>
        /// Returns the natural logarithm of the factorial of a number log(n!).
        /// </summary>
        /// <param name="n">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double LogFactorial(double n)
        {
            return Special.GammaLog(n + 1.0);
        }
        /// <summary>
        /// Returns the factorial of a number.
        /// </summary>
        /// <param name="n">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Factorial(double n)
        {
            // check it:
            if (Maths.IsInteger(n) && n >= 0)
            {
                // get it from memory
                if (n <= 170)
                {
                    return Special.A000142[(int)n];
                }
                return double.NaN;
            }
            return Special.Gamma(n + 1);
        }
        /// <summary>
        /// Returns the decreasing factorial of a number.
        /// </summary>
        /// <param name="n">Number</param>
        /// <param name="k">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double FactorialDown(double n, double k)
        {
            return Special.Factorial(n) / Special.Factorial(n - k);
        }
        /// <summary>
        /// Returns the increasing factorial of a number (Pohhammer symbol).
        /// </summary>
        /// <param name="n">Number</param>
        /// <param name="k">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double FactorialUp(double n, double k)
        {
            if (n == 0)
            {
                return 1.0;
            }
            return Special.Gamma(n + k) / Special.Gamma(n);
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
            double r = 2.2360679774997896964;
            double phi = (1.0 + r) / 2.0;
            double psi = (1.0 - r) / 2.0;
            double num = Math.Pow(phi, n) - Math.Pow(psi, n);
            return (int)(num / r);
        }
        /// <summary>
        /// Returns the value of the Luca number.
        /// </summary>
        /// <param name="n">Integer number</param>
        /// <returns>Integer number</returns>
        public static int Lucas(int n)
        {
            double r = 2.2360679774997896964;
            double phi = (1.0 + r) / 2.0;
            double psi = (1.0 - r) / 2.0;
            double num = Math.Pow(phi, n) + Math.Pow(psi, n);
            return (int)(num);
        }
        #endregion

        #region Bernoulli function
        /// <summary>
        /// Returns the Bernoulli number.
        /// </summary>
        /// <param name="n">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Bernoulli(int n)
        {
            // special cases:
            if (n < 0)
                return double.NaN;
            else if (n == 0)
                return 1;
            else if (n == 1)
                return -0.5;

            // for even number:
            else if (Maths.Mod(n, 2) == 0)
            {
                // get it from memory
                if (n <= 258)
                {
                    return Special.A027641[n / 2 - 1];
                }
                return double.NaN;
            }
            // for odd number:
            return 0;
        }
        /// <summary>
        /// Returns the value of the Bernoulli polynomial.
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Bernoulli(int n, double x)
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

        #region Euler function
        /// <summary>
        /// Returns the Euler number.
        /// </summary>
        /// <param name="n">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Euler(int n)
        {
            // special cases:
            if (n < 0)
                return double.NaN;
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
                return double.NaN;
            }
            // for odd number:
            return 0;
        }
        /// <summary>
        /// Returns the value of the Euler polynomial.
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="x">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Euler(int n, double x)
        {
            // properties:
            double p = 1, s = 0;
            double v = x - 0.5;
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

        #region Harmonic number
        /// <summary>
        /// Returns the harmonic number.
        /// </summary>
        /// <param name="n">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Harm(int n)
        {
            return Special.DiGamma(n) + 1.0 / n + Maths.Gamma;
        }
        /// <summary>
        /// Returns the harmonic number.
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="m">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Harm(int n, double m)
        {
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += Math.Pow(i + 1, -m);
            }

            return sum;
        }
        #endregion

        #region Chebyshev polynomials
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the first kind.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="n">Order</param>
        /// <returns>Double precision floating point number</returns>
        public static double ChebyshevT(double x, int n)
        {
            return Math.Cos(n * Math.Acos(x));
        }
        /// <summary>
        /// Returns the value of the Chebyshev polynomial of the second kind.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="n">Order</param>
        /// <returns>Double precision floating point number</returns>
        public static double ChebyshevU(double x, int n)
        {
            double z = Math.Acos(x);
            return Math.Sin((n + 1) * z) / Math.Sin(z);
        }
        #endregion

        #region Abel polynomials
        /// <summary>
        /// Returns the value of the Abel polynomial.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Power</param>
        /// <param name="n">Order</param>
        /// <returns>Double precision floating point number</returns>
        public static double Abel(double x, double a, int n)
        {
            if (n == 0)
                return 1;
            if (n == 1)
                return x;

            // Generalized formula
            // Abel polynomials recurrence relation for any n ≥ 1:
            return x * Math.Pow(x - a * n, n - 1);
        }
        #endregion

        #region Laguerre polynomial
        /// <summary>
        /// Returns the value of the Laguerre polynomial.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Power</param>
        /// <param name="k">Order</param>
        /// <returns>Double precision floating point number</returns>
        public static double Laguerre(double x, double a, int k)
        {
            if (k == 0)
                return 1;
            if (k == 1)
                return 1 + a - x;

            // Generalized formula
            // Laguerre polynomials recurrence relation for any k ≥ 1:
            double psi = (2 * k + 1 + a - x) * Laguerre(x, a, k - 1) - (k + a) * Laguerre(x, a, k - 2);
            double ksi = k + 1;
            return psi / ksi;
        }
        #endregion

        #region Legendre polynomial
        /// <summary>
        /// Returns the value of the Legendre polynomial of the first kind.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="m">Order</param>
        /// <returns>Double precision floating point number</returns>
        public static double Legendre(double x, int m)
        {
            if (m == 0)
                return 1;
            if (m == 1)
                return x;

            // legendre function
            // More info: https://pdfs.semanticscholar.org/be30/38b3aafed73d33fe78a701ee096471d16fa1.pdf
            // Laguerre polynomials recurrence relation for any k ≥ 1:
            double ksi = (2 * m + 1) * x * Legendre(x, m - 1) - m * Legendre(x, m - 2);
            double psi = m + 1;
            return ksi / psi;
        }
        #endregion

        #region Hermite polynomial
        /// <summary>
        /// Returns the value of the Hermite polynomial.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="m">Order</param>
        /// <returns>Double precision floating point number</returns>
        public static double Hermite(double x, int m)
        {
            if (m == 0)
                return 1;
            if (m == 1)
                return x;

            // recursion formula for Hermite polynomials
            double ksi = x * Hermite(x, m - 1) - m * Hermite(x, m - 2);
            return 2 * ksi;
        }
        #endregion

        #region Gegenbauer polynomial
        /// <summary>
        /// Returns the value of the Gegenbauer polynomial.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Power</param>
        /// <param name="n">Order</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gegenbauer(double x, double a, int n)
        {
            if (n == 0)
                return 1;
            if (n == 1)
                return 2 * a * x;

            // Generalized formula
            // Laguerre polynomials recurrence relation for any k ≥ 1:
            double psi = 2.0 * x * (n + a - 1) * Gegenbauer(x, a, n - 1) - (n + 2 * a - 2) * Gegenbauer(x, a, n - 2);
            return psi / n;
        }
        #endregion

        #region Mahlerpolynom function
        /// <summary>
        /// Returns the value of the Mahler polynomial.
        /// </summary>
        /// <param name="x">Number</param>
        /// <param name="t">Parameter</param>
        /// <returns>Double precision floating point number</returns>
        public static double Mahler(double x, double t)
        {
            return Math.Exp(x * (1 + t - Math.Pow(Math.E, t)));
        }
        #endregion

        #region Gompertz function
        /// <summary>
        /// Gets the value of the Gompertz function.
        /// </summary>
        /// <param name="t">Argument</param>
        /// <param name="a">Upper asymptote</param>
        /// <param name="b">Growth parameter</param>
        /// <param name="c">Growth rate</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gompertz(double t, double a, double b, double c)
        {
            double x = -c * t;
            double y = -b * Maths.E;
            double z = a * Maths.E;
            return Math.Pow(z, Math.Pow(y, x));
        }
        #endregion

        #region Heavyside delta-function
        /// <summary>
        /// Returns the value of the Heaviside delta function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="k">Smoothing factor</param>
        /// <returns>Double precision floating point number</returns>
        public static double Heaviside(double x, double k)
        {
            return 0.5 + 0.5 * Math.Tanh(k * x);
        }
        #endregion

        #region Logistic function
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Lower asymptote</param>
        /// <param name="k">Upper asymptote</param>
        /// <param name="b">Growth rate</param>
        /// <param name="v">Affect</param>
        /// <param name="q">Central moment</param>
        /// <param name="c">Offset</param>
        /// <returns>Double precision floating point number</returns>
        public static double Logistic(double x, double a, double k, double b, double v, double q, double c)
        {
            return a + (k - a) / Math.Pow(c + q * Math.Exp(-b * x), 1.0 / v);
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Lower asymptote</param>
        /// <param name="k">Upper asymptote</param>
        /// <param name="b">Growth rate</param>
        /// <returns>Double precision floating point number</returns>
        public static double Logistic(double x, double a, double k, double b)
        {
            return Logistic(x, a, k, b, 0.5, 0.5, 1);
        }
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Logistic(double x)
        {
            return Logistic(x, 0, 1, 3, 0.5, 0.5, 1);
        }
        #endregion

        #region Dirac delta-function
        /// <summary>
        /// Returns the value of the Dirac delta function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Coefficient</param>
        /// <returns>Double precision floating point number</returns>
        public static double Dirac(double x, double a)
        {
            double s = Math.Sqrt(Math.PI);
            double b = 1.0 / Math.Abs(a) / s;
            double c = Math.Pow(x / a, 2);
            double e = Math.Exp(-c);
            return b * e;
        }
        #endregion

        #region Dawson functions
        /// <summary>
        /// Returns the value of the D- / D + Dawson function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="positive">D- или D+</param>
        /// <returns>Double precision floating point number</returns>
        public static double Dawson(double x, bool positive)
        {
            if (positive)
            {
                // D+ function:
                double p = Special.sqrtPI / 2.0;
                double v = Maths.Exp(-x * x);
                double erfi = Special.Erfi(x);
                return p * v * erfi;
            }
            // D- function:
            double y = x * x;
            double g = sqrtPI / 2.0;
            double e = Special.Erf(x);
            double d = Math.Exp(y);
            return g * d * e;
        }
        #endregion

        #region Hypergeometric functions
        /// <summary>
        /// Returns the value of a hypergeometric function.
        /// <remarks>
        /// This version of the hypergeometric function is found in the Russian literature and is indicated: F(a,b,c,z).
        /// More information can be found on the website:
        /// https://en.wikipedia.org/wiki/Hypergeometric_function
        /// </remarks>
        /// </summary>
        /// <param name="a">Parameter</param>
        /// <param name="b">Parameter</param>
        /// <param name="c">Parameter</param>
        /// <param name="z">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Hypergeom(double a, double b, double c, double z)
        {
            // catch exception:
            if (double.IsNaN(a) || double.IsNaN(b) || double.IsNaN(c))
                throw new Exception("For this function, all input parameters must be correctly defined");

            // for all z = 0:
            if (z == 0)
                return 1;

            // Properties:
            double s = 1.0;
            double m = 1.0;
            double pa = 1, pb = 1, pc = 1;
            double t, eps = 1e-32;
            int i, j, iterations = 140;

            // Taylor series:
            for (i = 1; i < iterations; i++)
            {
                // Pochhammer symbols:
                j = i - 1;
                pa *= (a + j);
                pb *= (b + j);
                pc *= (c + j);

                // value:
                m *= z / i;
                t = (pa * pb * m) / pc;

                // stop point:
                if (Math.Abs(t) < eps)
                { break; }
                else { s += t; }
            }

            // result:
            return s;
        }
        /// <summary>
        /// Returns the value of a hypergeometric function.
        /// <remarks>
        /// The hypergeometric function can be used in several variations:
        /// F(a,b,z); F(a,~,z); F(~,b,z); F(~,~,z).
        /// Instead of the “~” sign, use the double.NaN value.
        /// More information can be found on the website:
        /// https://www.mathworks.com/help/symbolic/hypergeom.html#bt1nkmw-2
        /// </remarks>
        /// </summary>
        /// <param name="a">Parameter</param>
        /// <param name="b">Parameter</param>
        /// <param name="z">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Hypergeom(double a, double b, double z)
        {
            // for all z = 0:
            if (z == 0)
                return 1;

            // Properties:
            double s = 1.0;
            double m = 1.0;
            double pa = 1, pb = 1;
            double t, eps = 1e-32;
            int i, j, iterations = 140;

            if (double.IsNaN(a) && double.IsNaN(b) ||
                a == b)
            {
                // s = e^z:
                s = Math.Exp(z);
            }
            else if (double.IsNaN(b))
            {
                // special case for complex
                // values:
                if (Math.Abs(z) > 1.0)
                    return double.NaN;

                // F(a,~,z):
                for (i = 1; i < iterations; i++)
                {
                    // Pochhammer symbols:
                    j = i - 1;
                    pa *= (a + j);

                    // value:
                    m *= z / i;
                    t = (pa * m);

                    // stop point:
                    if (Math.Abs(t) < eps)
                    { break; }
                    else { s += t; }
                }
            }
            else if (double.IsNaN(a))
            {
                // F(~,b,z):
                for (i = 1; i < iterations; i++)
                {
                    // Pochhammer symbols:
                    j = i - 1;
                    pb *= (b + j);

                    // value:
                    m *= z / i;
                    t = m / pb;

                    // stop point:
                    if (Math.Abs(t) < eps)
                    { break; }
                    else { s += t; }
                }
            }
            else
            {
                // F(a,b,z):
                for (i = 1; i < iterations; i++)
                {
                    // Pochhammer symbols:
                    j = i - 1;
                    pa *= (a + j);
                    pb *= (b + j);

                    // value:
                    m *= z / i;
                    t = (pa * m) / pb;

                    // stop point:
                    if (Math.Abs(t) < eps)
                    { break; }
                    else { s += t; }
                }
            }

            // result:
            return s;
        }
        #endregion

        #region Generalized error function
        /// <summary>
        /// Returns the value of the generalized error function.
        /// </summary>
        /// <param name="x">Argument (0, +inf)</param>
        /// <param name="n">Order [0, +inf)</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gerf(double x, int n)
        {
            // Generalized error functions:
            if (n < 0)
                return double.NaN; // singular values

            else if (n == 0)
                return x / Math.E / sqrtPI; // E0(x)

            else if (n == 1)
                return (1 - Math.Exp(-x)) / sqrtPI; // E1(x)

            else if (n == 2)
                return Erf(x); // E2(x) = erf(x)

            // En(x) for all x > 0
            if (x > 0)
            {
                double p = 1.0 / n;
                double w = 1.0 / sqrtPI;
                double v = Gamma(p) - GammaIncomplete(p, Math.Pow(x, n), true);
                return w * Gamma(n) * (v);
            }
            return double.NaN;
        }
        /// <summary>
        /// Returns the value of the generalized error function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gerf(double x)
        {
            return Gerf(x, 2);
        }
        #endregion

        #region Rademacher function
        /// <summary>
        /// Returns the value of the Radamecher function.
        /// </summary>
        /// <param name="t">Argument [0, 1]</param>
        /// <param name="n">Order</param>
        /// <returns>Double precision floating point number</returns>
        public static double Rademacher(double t, int n)
        {
            double p = Math.Pow(2, n);
            double v = p * Math.PI * t;
            return Math.Sign(Math.Sin(v));
        }
        #endregion

        #region Elrang B and C functions
        /// <summary>
        /// Returns the value of the Erlang C-function.
        /// </summary>
        /// <param name="y">Firset parameter</param>
        /// <param name="v">Second parameter</param>
        /// <param name="t">Time parameter</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erlang(double y, int v, double t)
        {
            double e = Special.Erlang(y, v);
            double a = v * e;
            double b = v - y + y * e;
            double c = (v - y) * t;
            return a / b * Math.Exp(-c);
        }
        /// <summary>
        /// Returns the value of the Erlang B-function.
        /// </summary>
        /// <param name="y">Firset parameter</param>
        /// <param name="v">Second parameter</param>
        /// <returns>Double precision floating point number</returns>
        public static double Erlang(double y, int v)
        {
            // special cases:
            if (v == 0)
                return 1;
            if (v < 0)
                return double.NaN;

            // set:
            double t = 1, b = 1; int i;

            //series:
            for (i = 1; i < v; i++)
            {
                t *= y / i;
                b += t;
            }

            // last step and result:
            double a = (t * y / i);
            return a / (a + b);
        }
        #endregion

        #region Lambert W-function
        /// <summary>
        /// Returns the value of the Lambert W-function.
        /// </summary>
        /// <param name="x">Argument [-1/e,+inf)</param>
        /// <param name="branch">Function branch</param>
        /// <returns>Double precision floating point number</returns>
        public static double LambertW(double x, bool branch = true)
        {
            // *****************************************************************
            // Halley's method for calculating Lambert W-function.
            // *****************************************************************
            // Halley's method is a method for finding a zero of a real-valued 
            // function with a continuous and easily computed second derivative.

            double eps = 1e-8, e, f;
            double w = branch ? 1 : -2;
            double v = double.PositiveInfinity * w;

            while (Math.Abs(w - v) / Math.Abs(w) > eps)
            {
                v = w;
                e = Math.Exp(w);
                f = w * e - x;  // Iterate to make this quantity zero
                w = w - f / ((e * (w + 1) - (w + 2) * f / (2 * w + 2)));
            }
            return w;
        }
        /// <summary>
        /// Returns the value of the square super-root.
        /// </summary>
        /// <param name="x">Argument [1,+inf)</param>
        /// <param name="branch">Function branch</param>
        /// <returns>Double precision floating point number</returns>
        public static double Ssqrt(double x, bool branch = true)
        {
            // The 2nd-order super-root, square super-root, or super square root has notation ssqrt(x).
            // It can be represented with the Lambert W-function: ssqrt(x) = log(x) / W{ log(x) }.
            double log = Math.Log(x);
            return log / LambertW(log, branch);
        }
        #endregion

        #region Minkowski function
        /// <summary>
        /// Returns the value of the Minkowski function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Minkowski(double x)
        {
            // Minkowski function:
            long p = (long)x, q = 1, r = p + 1, s = 1, m, n;
            double d = 1, y = p;

            if (x < (double)p || (p < 0) ^ (r <= 0))
            {
                return x;
            }

            // calculating:
            for (; ; )
            {
                d /= 2; if (y + d == y) break;
                m = p + r; if ((m < 0) ^ (p < 0)) break;
                n = q + s; if (n < 0) break;

                if (x < (double)m / n)
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
        /// 
        /// </summary>
        private static double[] A027641 = new double[]
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
        /// 
        /// </summary>
        private static double[] A000142 = new double[]
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
        /// 
        /// </summary>
        private static double[] A122045 = new double[]
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
