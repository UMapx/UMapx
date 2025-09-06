using System;
using System.Text;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements the least squares method.
    /// </summary>
    internal static class LeastSquaresOptions
    {
        #region Float components
        /// <summary>
        /// Returns the polynomial value.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="c">Approximation coefficients</param>
        /// <returns>Value</returns>
        public static float Polynomial(float x, float[] c)
        {
            int n = c.Length, i;
            float p = 1, s = 0;

            for (i = 0; i < n; i++, p *= x)
            {
                s += c[i] * p;
            }
            return s;
        }
        /// <summary>
        /// Returns an array of polynomial values.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="c">Approximation coefficients</param>
        /// <returns>Array</returns>
        public static float[] Polynomial(float[] x, float[] c)
        {
            int n = x.Length, i;
            float[] y = new float[n];

            for (i = 0; i < n; i++)
            {
                y[i] = LeastSquaresOptions.Polynomial(x[i], c);
            }
            return y;
        }
        /// <summary>
        /// Returns an array of polynomial values.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="y">Function</param>
        /// <param name="iterations">Number of iterations</param>
        /// <returns>Array</returns>
        public static float[] Coefficients(float[] x, float[] y, int iterations)
        {
            int i, j;
            int m = iterations < 1 ? 1 : iterations;
            float[,] matrix = new float[m, m + 1];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    matrix[i, j] = LeastSquaresOptions.SummaryPow(x, j + i);
                }
                matrix[i, m] = LeastSquaresOptions.SummaryPow(y, x, 1, i);
            }

            return Matrice.Solve(matrix);
        }
        /// <summary>
        /// Returns the value of the expression: s += v(i) ^ pow.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="pow">Power</param>
        /// <returns>Value</returns>
        public static float SummaryPow(float[] v, float pow)
        {
            float sum = 0;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                sum += Maths.Pow(v[i], pow);
            }
            return sum;
        }
        /// <summary>
        /// Returns the value of the expression: s += {x(i) ^ powx} * {y(i) ^ powy}.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="y">Array</param>
        /// <param name="powx">Power of x</param>
        /// <param name="powy">Power of y</param>
        /// <returns>Value</returns>
        public static float SummaryPow(float[] x, float[] y, float powx, float powy)
        {
            float sum = 0;
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                sum += Maths.Pow(x[i], powx) * Maths.Pow(y[i], powy);
            }
            return sum;
        }
        /// <summary>
        /// Returns the approximation error of the function.
        /// </summary>
        /// <param name="a">Approximation</param>
        /// <param name="b">Function</param>
        /// <returns>Value</returns>
        public static float Error(float[] a, float[] b)
        {
            float vara = Matrice.Var(a);
            float varb = Matrice.Var(b);

            if (vara < varb)
            {
                return vara / varb;
            }
            return varb / vara;
        }
        /// <summary>
        /// Returns the equation of a polynomial represented as a string (uses " * X^" for powers).
        /// </summary>
        /// <param name="p">Polynomial coefficients</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Equation(float[] p) => Equation(p, " * X^");

        /// <summary>
        /// Returns the equation of a polynomial represented as a string.
        /// The <paramref name="function"/> string is appended for i>0 before the power index (e.g. " * X^").
        /// </summary>
        /// <param name="p">Polynomial coefficients</param>
        /// <param name="function">Token placed before the power index for i&gt;0 (e.g. " * X^")</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Equation(float[] p, string function)
        {
            if (p == null) throw new ArgumentNullException(nameof(p));
            if (function == null) throw new ArgumentNullException(nameof(function));

            int n = p.Length;
            if (n == 0) return string.Empty;

            // Preallocate roughly: coeff + token + digits + separators
            var sb = new StringBuilder(n * (function.Length + 8));

            for (int i = 0; i < n; i++)
            {
                sb.Append(p[i].ToString());                // keep current-culture formatting (same as original)
                if (i > 0) sb.Append(function).Append(i);  // append "function + exponent" for i>0

                if (i < n - 1)
                    sb.Append(p[i + 1] < 0f ? " " : " + "); // same sign/spacing rule as original
            }

            return sb.ToString();
        }
        #endregion

        #region Complex components
        /// <summary>
        /// Returns the polynomial value.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="c">Approximation coefficients</param>
        /// <returns>Complex number</returns>
        public static Complex32 Polynomial(Complex32 x, Complex32[] c)
        {
            int n = c.Length, i;
            Complex32 p = 1, s = 0;

            for (i = 0; i < n; i++, p *= x)
            {
                s += c[i] * p;
            }
            return s;
        }
        /// <summary>
        /// Returns an array of polynomial values.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="c">Approximation coefficients</param>
        /// <returns>Array</returns>
        public static Complex32[] Polynomial(Complex32[] x, Complex32[] c)
        {
            int n = x.Length, i;
            Complex32[] y = new Complex32[n];

            for (i = 0; i < n; i++)
            {
                y[i] = LeastSquaresOptions.Polynomial(x[i], c);
            }
            return y;
        }
        /// <summary>
        /// Returns an array of polynomial values.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="y">Function</param>
        /// <param name="iterations">Number of iterations</param>
        /// <returns>Array</returns>
        public static Complex32[] Coefficients(Complex32[] x, Complex32[] y, int iterations)
        {
            int i, j;
            int m = iterations < 1 ? 1 : iterations;
            Complex32[,] matrix = new Complex32[m, m + 1];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    matrix[i, j] = LeastSquaresOptions.SummaryPow(x, j + i);
                }
                matrix[i, m] = LeastSquaresOptions.SummaryPow(y, x, 1, i);
            }

            return Matrice.Solve(matrix);
        }
        /// <summary>
        /// Returns the value of the expression: s += v(i) ^ pow.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="pow">Power</param>
        /// <returns>Value</returns>
        public static Complex32 SummaryPow(Complex32[] v, float pow)
        {
            Complex32 sum = 0;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                sum += Maths.Pow(v[i], pow);
            }
            return sum;
        }
        /// <summary>
        /// Returns the value of the expression: s += {x(i) ^ powx} * {y(i) ^ powy}.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="y">Array</param>
        /// <param name="powx">Power of x</param>
        /// <param name="powy">Power of y</param>
        /// <returns>Value</returns>
        public static Complex32 SummaryPow(Complex32[] x, Complex32[] y, float powx, float powy)
        {
            Complex32 sum = 0;
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                sum += Maths.Pow(x[i], powx) * Maths.Pow(y[i], powy);
            }
            return sum;
        }
        /// <summary>
        /// Returns the approximation error of the function.
        /// </summary>
        /// <param name="a">Approximation</param>
        /// <param name="b">Function</param>
        /// <returns>Value</returns>
        public static Complex32 Error(Complex32[] a, Complex32[] b)
        {
            Complex32 vara = Matrice.Var(a);
            Complex32 varb = Matrice.Var(b);

            if (vara.Abs < varb.Abs)
            {
                return (vara / varb).Real;
            }
            return (varb / vara).Real;
        }
        /// <summary>
        /// Returns the equation of a polynomial represented as a string (uses " * X^" for powers).
        /// </summary>
        /// <param name="p">Polynomial coefficients</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Equation(Complex32[] p) => Equation(p, " * X^");

        /// <summary>
        /// Returns the equation of a polynomial represented as a string.
        /// The <paramref name="function"/> string is appended for i>0 before the power index (e.g. " * X^").
        /// </summary>
        /// <param name="p">Polynomial coefficients</param>
        /// <param name="function">Token placed before the power index for i&gt;0 (e.g. " * X^")</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Equation(Complex32[] p, string function)
        {
            if (p == null) throw new ArgumentNullException(nameof(p));
            if (function == null) throw new ArgumentNullException(nameof(function));

            int n = p.Length;
            if (n == 0) return string.Empty;

            // Preallocate roughly: coeff + token + digits + separators
            var sb = new StringBuilder(n * (function.Length + 8));

            for (int i = 0; i < n; i++)
            {
                sb.Append(p[i].ToString());                // keep current-culture formatting (same as original)
                if (i > 0) sb.Append(function).Append(i);  // append "function + exponent" for i>0

                if (i < n - 1)
                    sb.Append(p[i + 1].Abs < 0f ? " " : " + "); // same sign/spacing rule as original
            }

            return sb.ToString();
        }
        #endregion
    }
}
