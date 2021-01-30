using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines the least squares approximation class.
    /// <remarks>
    /// This class is a solution to the problem of finding the function A (x) ≈ F (x), where F (x) is the original function.
    /// More information can be found on the website:
    /// http://simenergy.ru/math-analysis/digital-processing/85-ordinary_least_squares
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Approximation
    {
        #region Private data
        private Approximation.Method method;
        private int power;
        #endregion

        #region Approximation components
        /// <summary>
        /// Initializes the least squares approximation class.
        /// </summary>
        /// <param name="power">Polynomial degree</param>
        /// <param name="method">Approximation method</param>
        public Approximation(int power = 1, Approximation.Method method = Approximation.Method.Polynomial)
        {
            this.Power = power;
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the degree of the polynomial.
        /// </summary>
        public int Power
        {
            get
            {
                return this.power;
            }
            set
            {
                if (value < 1)
                    throw new Exception("Invalid argument value");

                this.power = value;
            }
        }
        /// <summary>
        /// Gets or sets the approximation method.
        /// </summary>
        public Approximation.Method MethodType
        {
            get
            {
                return this.method;
            }
            set
            {
                this.method = value;
            }
        }
        #endregion

        #region Public voids
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <returns>Array</returns>
        public double[] Compute(double[] x, double[] y)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, out _, out _, out _);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, out _, out _, out _);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, out _, out _, out _);

                default:
                    return Approximation.powr(x, y, power, out _, out _, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <returns>Array</returns>
        public double[] Compute(double[] x, double[] y, out double[] cf)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out _, out _);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out _, out _);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, out cf, out _, out _);

                default:
                    return Approximation.powr(x, y, power, out cf, out _, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="similarity">Similarity</param>
        /// <returns>Array</returns>
        public double[] Compute(double[] x, double[] y, out double[] cf, out double similarity)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out similarity, out _);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out similarity, out _);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, out cf, out similarity, out _);

                default:
                    return Approximation.powr(x, y, power, out cf, out similarity, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="similarity">Similarity</param>
        /// <param name="equation">Equation</param>
        /// <returns>Array</returns>
        public double[] Compute(double[] x, double[] y, out double[] cf, out double similarity, out string equation)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out similarity, out equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out similarity, out equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, out cf, out similarity, out equation);

                default:
                    return Approximation.powr(x, y, power, out cf, out similarity, out equation);
            }
        }

        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <returns>Array</returns>
        public Complex[] Compute(Complex[] x, Complex[] y)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, out _, out _, out _);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, out _, out _, out _);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, out _, out _, out _);

                default:
                    return Approximation.powr(x, y, power, out _, out _, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <returns>Array</returns>
        public Complex[] Compute(Complex[] x, Complex[] y, out Complex[] cf)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out _, out _);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out _, out _);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, out cf, out _, out _);

                default:
                    return Approximation.powr(x, y, power, out cf, out _, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="similarity">Similarity</param>
        /// <returns>Array</returns>
        public Complex[] Compute(Complex[] x, Complex[] y, out Complex[] cf, out Complex similarity)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out similarity, out _);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out similarity, out _);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, out cf, out similarity, out _);

                default:
                    return Approximation.powr(x, y, power, out cf, out similarity, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="similarity">Similarity</param>
        /// <param name="equation">Equation</param>
        /// <returns>Array</returns>
        public Complex[] Compute(Complex[] x, Complex[] y, out Complex[] cf, out Complex similarity, out string equation)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out similarity, out equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out similarity, out equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, out cf, out similarity, out equation);

                default:
                    return Approximation.powr(x, y, power, out cf, out similarity, out equation);
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="power"></param>
        /// <param name="cf"></param>
        /// <param name="error"></param>
        /// <param name="equation"></param>
        /// <returns></returns>
        private static double[] poly(double[] x, double[] y, int power, out double[] cf, out double error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            cf = LeastSquaresOptions.Coefficients(x, y, m);
            double[] ya = LeastSquaresOptions.Polynomial(x, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf);
            return ya;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="power"></param>
        /// <param name="cf"></param>
        /// <param name="error"></param>
        /// <param name="equation"></param>
        /// <returns></returns>
        private static Complex[] poly(Complex[] x, Complex[] y, int power, out Complex[] cf, out Complex error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            cf = LeastSquaresOptions.Coefficients(x, y, m);
            Complex[] ya = LeastSquaresOptions.Polynomial(x, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf);
            return ya;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="power"></param>
        /// <param name="cf"></param>
        /// <param name="error"></param>
        /// <param name="equation"></param>
        /// <returns></returns>
        private static double[] logc(double[] x, double[] y, int power, out double[] cf, out double error, out string equation)
        {
            // Options:
            int n = x.Length, i;
            int m = (power < 1) ? 2 : power + 1;
            double[] xa = new double[n];
            double[] ya = new double[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Maths.Log(x[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, y, m);
            ya = LeastSquaresOptions.Polynomial(xa, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf, " * LN(X)^");
            return ya;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="power"></param>
        /// <param name="cf"></param>
        /// <param name="error"></param>
        /// <param name="equation"></param>
        /// <returns></returns>
        private static Complex[] logc(Complex[] x, Complex[] y, int power, out Complex[] cf, out Complex error, out string equation)
        {
            // Options:
            int n = x.Length, i;
            int m = (power < 1) ? 2 : power + 1;
            Complex[] xa = new Complex[n];
            Complex[] ya = new Complex[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Maths.Log(x[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, y, m);
            ya = LeastSquaresOptions.Polynomial(xa, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf, " * LN(X)^");
            return ya;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="power"></param>
        /// <param name="cf"></param>
        /// <param name="error"></param>
        /// <param name="equation"></param>
        /// <returns></returns>
        private static double[] expn(double[] x, double[] y, int power, out double[] cf, out double error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            double[] ya = new double[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Log(y[i], Math.E);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(x, ya, m);
            double[] p = LeastSquaresOptions.Polynomial(x, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Pow(Math.E, p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = "EXP" + '(' + LeastSquaresOptions.Equation(cf) + ')';
            return ya;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="power"></param>
        /// <param name="cf"></param>
        /// <param name="error"></param>
        /// <param name="equation"></param>
        /// <returns></returns>
        private static Complex[] expn(Complex[] x, Complex[] y, int power, out Complex[] cf, out Complex error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            Complex[] ya = new Complex[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Log(y[i], Math.E);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(x, ya, m);
            Complex[] p = LeastSquaresOptions.Polynomial(x, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Pow(Math.E, p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = "EXP" + '(' + LeastSquaresOptions.Equation(cf) + ')';
            return ya;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="power"></param>
        /// <param name="cf"></param>
        /// <param name="error"></param>
        /// <param name="equation"></param>
        /// <returns></returns>
        private static double[] powr(double[] x, double[] y, int power, out double[] cf, out double error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            double[] xa = new double[n];
            double[] ya = new double[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Maths.Log(x[i]);
                ya[i] = Maths.Log(y[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, ya, m);
            double[] p = LeastSquaresOptions.Polynomial(xa, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Exp(p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = "EXP" + '(' + LeastSquaresOptions.Equation(cf, " * LN(X)^") + ')';
            return ya;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="power"></param>
        /// <param name="cf"></param>
        /// <param name="error"></param>
        /// <param name="equation"></param>
        /// <returns></returns>
        private static Complex[] powr(Complex[] x, Complex[] y, int power, out Complex[] cf, out Complex error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            Complex[] xa = new Complex[n];
            Complex[] ya = new Complex[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Maths.Log(x[i]);
                ya[i] = Maths.Log(y[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, ya, m);
            Complex[] p = LeastSquaresOptions.Polynomial(xa, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Exp(p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = "EXP" + '(' + LeastSquaresOptions.Equation(cf, " * LN(X)^") + ')';
            return ya;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Approximation method.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Polynomial approximation.
            /// </summary>
            Polynomial,
            /// <summary>
            /// Logarithmic approximation.
            /// </summary>
            Logarithmic,
            /// <summary>
            /// Exponential approximation.
            /// </summary>
            Exponential,
            /// <summary>
            /// Power approximation.
            /// </summary>
            Power,
            #endregion
        }
        #endregion
    }
}
