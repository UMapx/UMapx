using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements the least squares method.
    /// </summary>
    internal static class LeastSquaresOptions
    {
        #region double components
        /// <summary>
        /// Returns the polynomial value.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="c">Approximation coefficients</param>
        /// <returns>Double precision floating point number</returns>
        public static double Polynomial(double x, double[] c)
        {
            int n = c.Length, i;
            double p = 1, s = 0;

            for (i = 0; i < n; i++, p *= x)
            {
                s += c[i] * p;
            }
            return s;
        }
        /// <summary>
        /// Returns an array of polynomial values.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="c">Approximation coefficients</param>
        /// <returns>Array</returns>
        public static double[] Polynomial(double[] x, double[] c)
        {
            int n = x.Length, i;
            double[] y = new double[n];

            for (i = 0; i < n; i++)
            {
                y[i] = LeastSquaresOptions.Polynomial(x[i], c);
            }
            return y;
        }
        /// <summary>
        /// Returns an array of polynomial values.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="y">Function</param>
        /// <param name="iterations">Number of iterations</param>
        /// <returns>Array</returns>
        public static double[] Coefficients(double[] x, double[] y, int iterations)
        {
            int i, j;
            int n = x.Length;
            int m = iterations < 1 ? 1 : iterations;
            double[,] matrix = new double[m, m + 1];

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
        /// <returns>Double precision floating point number</returns>
        public static double SummaryPow(double[] v, double pow)
        {
            double sum = 0;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                sum += Math.Pow(v[i], pow);
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
        /// <returns>Double precision floating point number</returns>
        public static double SummaryPow(double[] x, double[] y, double powx, double powy)
        {
            double sum = 0;
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                sum += Math.Pow(x[i], powx) * Math.Pow(y[i], powy);
            }
            return sum;
        }
        /// <summary>
        /// Returns the approximation error of the function.
        /// </summary>
        /// <param name="a">Approximation</param>
        /// <param name="b">Function</param>
        /// <returns>Double precision floating point number</returns>
        public static double Error(double[] a, double[] b)
        {
            double vara = Matrice.Var(a);
            double varb = Matrice.Var(b);

            if (vara < varb)
            {
                return vara / varb;
            }
            return varb / vara;
        }
        /// <summary>
        /// Returns the equation of a polynomial represented as a string.
        /// </summary>
        /// <param name="p">Polynomial coefficients</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Equation(double[] p)
        {
            string equation = "";
            int length = p.Length;

            for (int i = 0; i < length; i++)
            {
                equation += (Convert.ToString(p[i]) +
                            (i == 0 ? "" : (" * X^" + Convert.ToString(i))) +
                            (i < length - 1 ? (p[i + 1] < 0 ? " " : " + ") : ""));
            }

            return equation;
        }
        /// <summary>
        /// Returns the equation of a polynomial represented as a string.
        /// </summary>
        /// <param name="p">Polynomial coefficients</param>
        /// <param name="function">Function</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Equation(double[] p, string function)
        {
            string equation = "";
            int length = p.Length;

            for (int i = 0; i < length; i++)
            {
                equation += (Convert.ToString(p[i]) +
                            (i == 0 ? "" : (function + Convert.ToString(i))) +
                            (i < length - 1 ? (p[i + 1] < 0 ? " " : " + ") : ""));
            }

            return equation;
        }
        #endregion

        #region complex components
        /// <summary>
        /// Returns the polynomial value.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="c">Approximation coefficients</param>
        /// <returns>Complex number</returns>
        public static Complex Polynomial(Complex x, Complex[] c)
        {
            int n = c.Length, i;
            Complex p = 1, s = 0;

            for (i = 0; i < n; i++, p *= x)
            {
                s += c[i] * p;
            }
            return s;
        }
        /// <summary>
        /// Returns an array of polynomial values.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="c">Approximation coefficients</param>
        /// <returns>Array</returns>
        public static Complex[] Polynomial(Complex[] x, Complex[] c)
        {
            int n = x.Length, i;
            Complex[] y = new Complex[n];

            for (i = 0; i < n; i++)
            {
                y[i] = LeastSquaresOptions.Polynomial(x[i], c);
            }
            return y;
        }
        /// <summary>
        /// Returns an array of polynomial values.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="y">Function</param>
        /// <param name="iterations">Number of iterations</param>
        /// <returns>Array</returns>
        public static Complex[] Coefficients(Complex[] x, Complex[] y, int iterations)
        {
            int i, j;
            int n = x.Length;
            int m = iterations < 1 ? 1 : iterations;
            Complex[,] matrix = new Complex[m, m + 1];

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
        /// <returns>Double precision floating point number</returns>
        public static Complex SummaryPow(Complex[] v, double pow)
        {
            Complex sum = 0;
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
        /// <returns>Double precision floating point number</returns>
        public static Complex SummaryPow(Complex[] x, Complex[] y, double powx, double powy)
        {
            Complex sum = 0;
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
        /// <returns>Double precision floating point number</returns>
        public static Complex Error(Complex[] a, Complex[] b)
        {
            Complex vara = Matrice.Var(a);
            Complex varb = Matrice.Var(b);

            if (vara.Abs < varb.Abs)
            {
                return (vara / varb).Real;
            }
            return (varb / vara).Real;
        }
        /// <summary>
        /// Returns the equation of a polynomial represented as a string.
        /// </summary>
        /// <param name="p">Polynomial coefficients</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Equation(Complex[] p)
        {
            string equation = "";
            int length = p.Length;

            for (int i = 0; i < length; i++)
            {
                equation += ("(" + Convert.ToString(p[i]) + ")" +
                            (i == 0 ? "" : (" * X^" + Convert.ToString(i))) +
                            (i < length - 1 ? (p[i + 1].Abs < 0 ? " " : " + ") : ""));
            }

            return equation;
        }
        /// <summary>
        /// Returns the equation of a polynomial represented as a string.
        /// </summary>
        /// <param name="p">Polynomial coefficients</param>
        /// <param name="function">Function</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Equation(Complex[] p, string function)
        {
            string equation = "";
            int length = p.Length;

            for (int i = 0; i < length; i++)
            {
                equation += ("(" + Convert.ToString(p[i]) + ")" +
                            (i == 0 ? "" : (function + Convert.ToString(i))) +
                            (i < length - 1 ? (p[i + 1].Abs < 0 ? " " : " + ") : ""));
            }

            return equation;
        }
        #endregion
    }
}
