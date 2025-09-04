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

            return MatrixF.Solve(matrix);
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
                sum += (float)Math.Pow(v[i], pow);
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
                sum += (float)Math.Pow(x[i], powx) * (float)Math.Pow(y[i], powy);
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
            float vara = MatrixF.Var(a);
            float varb = MatrixF.Var(b);

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
        public static ComplexF Polynomial(ComplexF x, ComplexF[] c)
        {
            int n = c.Length, i;
            ComplexF p = 1, s = 0;

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
        public static ComplexF[] Polynomial(ComplexF[] x, ComplexF[] c)
        {
            int n = x.Length, i;
            ComplexF[] y = new ComplexF[n];

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
        public static ComplexF[] Coefficients(ComplexF[] x, ComplexF[] y, int iterations)
        {
            int i, j;
            int m = iterations < 1 ? 1 : iterations;
            ComplexF[,] matrix = new ComplexF[m, m + 1];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    matrix[i, j] = LeastSquaresOptions.SummaryPow(x, j + i);
                }
                matrix[i, m] = LeastSquaresOptions.SummaryPow(y, x, 1, i);
            }

            return MatrixF.Solve(matrix);
        }
        /// <summary>
        /// Returns the value of the expression: s += v(i) ^ pow.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="pow">Power</param>
        /// <returns>Value</returns>
        public static ComplexF SummaryPow(ComplexF[] v, float pow)
        {
            ComplexF sum = 0;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                sum += MathsF.Pow(v[i], pow);
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
        public static ComplexF SummaryPow(ComplexF[] x, ComplexF[] y, float powx, float powy)
        {
            ComplexF sum = 0;
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                sum += MathsF.Pow(x[i], powx) * MathsF.Pow(y[i], powy);
            }
            return sum;
        }
        /// <summary>
        /// Returns the approximation error of the function.
        /// </summary>
        /// <param name="a">Approximation</param>
        /// <param name="b">Function</param>
        /// <returns>Value</returns>
        public static ComplexF Error(ComplexF[] a, ComplexF[] b)
        {
            ComplexF vara = MatrixF.Var(a);
            ComplexF varb = MatrixF.Var(b);

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
        public static string Equation(ComplexF[] p) => Equation(p, " * X^");

        /// <summary>
        /// Returns the equation of a polynomial represented as a string.
        /// The <paramref name="function"/> string is appended for i>0 before the power index (e.g. " * X^").
        /// </summary>
        /// <param name="p">Polynomial coefficients</param>
        /// <param name="function">Token placed before the power index for i&gt;0 (e.g. " * X^")</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Equation(ComplexF[] p, string function)
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
