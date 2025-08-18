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
        private ApproximationMethod method;
        private int power;
        #endregion

        #region Approximation components
        /// <summary>
        /// Initializes the least squares approximation class.
        /// </summary>
        /// <param name="power">Polynomial degree</param>
        /// <param name="method">Approximation method</param>
        public Approximation(int power = 1, ApproximationMethod method = ApproximationMethod.Polynomial)
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
                    throw new ArgumentException("Invalid argument value");

                this.power = value;
            }
        }
        /// <summary>
        /// Gets or sets the approximation method.
        /// </summary>
        public ApproximationMethod MethodType
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
        public float[] Compute(float[] x, float[] y)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.poly(x, y, power, out _, out _, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.logc(x, y, power, out _, out _, out _);

                case ApproximationMethod.Exponential:
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
        public float[] Compute(float[] x, float[] y, out float[] cf)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out _, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out _, out _);

                case ApproximationMethod.Exponential:
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
        public float[] Compute(float[] x, float[] y, out float[] cf, out float similarity)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out similarity, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out similarity, out _);

                case ApproximationMethod.Exponential:
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
        public float[] Compute(float[] x, float[] y, out float[] cf, out float similarity, out string equation)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out similarity, out equation);

                case ApproximationMethod.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out similarity, out equation);

                case ApproximationMethod.Exponential:
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
        public Complex32[] Compute(Complex32[] x, Complex32[] y)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.poly(x, y, power, out _, out _, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.logc(x, y, power, out _, out _, out _);

                case ApproximationMethod.Exponential:
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
        public Complex32[] Compute(Complex32[] x, Complex32[] y, out Complex32[] cf)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out _, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out _, out _);

                case ApproximationMethod.Exponential:
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
        public Complex32[] Compute(Complex32[] x, Complex32[] y, out Complex32[] cf, out Complex32 similarity)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out similarity, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out similarity, out _);

                case ApproximationMethod.Exponential:
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
        public Complex32[] Compute(Complex32[] x, Complex32[] y, out Complex32[] cf, out Complex32 similarity, out string equation)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.poly(x, y, power, out cf, out similarity, out equation);

                case ApproximationMethod.Logarithmic:
                    return Approximation.logc(x, y, power, out cf, out similarity, out equation);

                case ApproximationMethod.Exponential:
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
        private static float[] poly(float[] x, float[] y, int power, out float[] cf, out float error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            cf = LeastSquaresOptions.Coefficients(x, y, m);
            float[] ya = LeastSquaresOptions.Polynomial(x, cf);
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
        private static Complex32[] poly(Complex32[] x, Complex32[] y, int power, out Complex32[] cf, out Complex32 error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            cf = LeastSquaresOptions.Coefficients(x, y, m);
            Complex32[] ya = LeastSquaresOptions.Polynomial(x, cf);
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
        private static float[] logc(float[] x, float[] y, int power, out float[] cf, out float error, out string equation)
        {
            // Options:
            int n = x.Length, i;
            int m = (power < 1) ? 2 : power + 1;
            float[] xa = new float[n];
            float[] ya;

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Maths.Log(x[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, y, m);
            ya = LeastSquaresOptions.Polynomial(xa, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf, " * Log(x)^");
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
        private static Complex32[] logc(Complex32[] x, Complex32[] y, int power, out Complex32[] cf, out Complex32 error, out string equation)
        {
            // Options:
            int n = x.Length, i;
            int m = (power < 1) ? 2 : power + 1;
            Complex32[] xa = new Complex32[n];
            Complex32[] ya;

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Maths.Log(x[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, y, m);
            ya = LeastSquaresOptions.Polynomial(xa, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf, " * Log(x)^");
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
        private static float[] expn(float[] x, float[] y, int power, out float[] cf, out float error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            float[] ya = new float[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Log(y[i], Maths.E);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(x, ya, m);
            float[] p = LeastSquaresOptions.Polynomial(x, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Pow(Maths.E, p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = $"Exp({LeastSquaresOptions.Equation(cf)})";
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
        private static Complex32[] expn(Complex32[] x, Complex32[] y, int power, out Complex32[] cf, out Complex32 error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            Complex32[] ya = new Complex32[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Log(y[i], Maths.E);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(x, ya, m);
            Complex32[] p = LeastSquaresOptions.Polynomial(x, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Pow(Maths.E, p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = $"Exp({LeastSquaresOptions.Equation(cf)})";
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
        private static float[] powr(float[] x, float[] y, int power, out float[] cf, out float error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            float[] xa = new float[n];
            float[] ya = new float[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Maths.Log(x[i]);
                ya[i] = Maths.Log(y[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, ya, m);
            float[] p = LeastSquaresOptions.Polynomial(xa, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Exp(p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = $"Exp({LeastSquaresOptions.Equation(cf, " * Log(x)^")})";
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
        private static Complex32[] powr(Complex32[] x, Complex32[] y, int power, out Complex32[] cf, out Complex32 error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            Complex32[] xa = new Complex32[n];
            Complex32[] ya = new Complex32[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Maths.Log(x[i]);
                ya[i] = Maths.Log(y[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, ya, m);
            Complex32[] p = LeastSquaresOptions.Polynomial(xa, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Exp(p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = $"Exp({LeastSquaresOptions.Equation(cf, " * Log(x)^")})";
            return ya;
        }
        #endregion
    }
}
