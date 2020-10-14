using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements the solution of a nonlinear equation.
    /// <remarks>
    /// This class is a solution to the problem of finding the root of a nonlinear equation of the form F(x) = 0.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Nonlinear
    {
        #region Private data
        private Nonlinear.Method method;
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements the solution of a nonlinear equation.
        /// </summary>
        /// <param name="eps">Epsilon [0, 1]</param>
        /// <param name="method">Method for solving a nonlinear equation</param>
        public Nonlinear(double eps = 1e-8, Nonlinear.Method method = Method.Secant)
        {
            this.method = method;
            this.Eps = eps;
        }
        /// <summary>
        /// Gets or sets the method for solving the nonlinear equation.
        /// </summary>
        public Nonlinear.Method MethodType
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
        /// <summary>
        /// Gets or sets the error value [0, 1].
        /// </summary>
        public double Eps
        {
            get
            {
                return this.eps;
            }
            set
            {
                this.eps = Maths.Double(value);
            }
        }
        /// <summary>
        /// Gets the root value of a nonlinear equation.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Start of line</param>
        /// <param name="b">End of line</param>
        /// <returns>Double precision floating point number</returns>
        public double Compute(IDouble function, double a, double b)
        {
            // chose method of nonlinear
            switch (method)
            {
                case Method.Chord:
                    return Nonlinear.chord(function, a, b, this.eps);
                case Method.FalsePosition:
                    return Nonlinear.falpo(function, a, b, this.eps);
                case Method.Secant:
                    return Nonlinear.secan(function, a, b, this.eps);

                default:
                    return Nonlinear.bisec(function, a, b, this.eps);
            }
        }
        /// <summary>
        /// Gets the root value of a nonlinear equation.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Start of line</param>
        /// <param name="b">End of line</param>
        /// <returns>Double precision floating point number</returns>
        public Complex Compute(IComplex function, Complex a, Complex b)
        {
            // chose method of nonlinear
            switch (method)
            {
                case Method.Chord:
                    return Nonlinear.chord(function, a, b, this.eps);
                case Method.FalsePosition:
                    return Nonlinear.falpo(function, a, b, this.eps);

                default:
                    return Nonlinear.secan(function, a, b, this.eps);
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static double bisec(IDouble f, double a, double b, double eps = 1e-8)
        {
            double x1 = a; double x2 = b;
            double fb = f(b);
            double midpt;
            int n = 0;

            while (Math.Abs(x2 - x1) > eps && n < short.MaxValue)
            {
                midpt = 0.5 * (x1 + x2);

                if (fb * f(midpt) > 0)
                    x2 = midpt;
                else
                    x1 = midpt;
                n++;
            }
            return x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static double secan(IDouble f, double a, double b, double eps = 1e-8)
        {
            double x1 = a;
            double x2 = b;
            double fb = f(b);
            double mpoint;
            int n = 0;

            while (Math.Abs(f(x2)) > eps && n < short.MaxValue)
            {
                mpoint = x2 - (x2 - x1) * fb / (fb - f(x1));
                x1 = x2;
                x2 = mpoint;
                fb = f(x2);
                n++;
            }
            return x2;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static double falpo(IDouble f, double a, double b, double eps = 1e-8)
        {
            double x1 = a;
            double x2 = b;
            double fb = f(b);
            int n = 0;

            while (Math.Abs(x2 - x1) > eps && n < short.MaxValue)
            {
                double xpoint = x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
                if (fb * f(xpoint) > 0)
                    x2 = xpoint;
                else
                    x1 = xpoint;
                if (Math.Abs(f(xpoint)) < eps)
                    break;
                n++;
            }
            return x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static double chord(IDouble f, double a, double b, double eps = 1e-8)
        {
            int n = 0;
            double x0 = (b - a) / 2.0;
            double x;

            while (Math.Abs(f(x0) / b) > eps && n < short.MaxValue)
            {
                x = x0;
                x0 = x - (f(x) * (a - x)) / (f(a) - f(x));
                n++;
            }
            return x0;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static Complex chord(IComplex f, Complex a, Complex b, double eps = 1e-8)
        {
            int n = 0;
            Complex x0 = (b - a) / 2.0;
            Complex x;

            while (Maths.Abs(f(x0) / b) > eps && n < short.MaxValue)
            {
                x = x0;
                x0 = x - (f(x) * (a - x)) / (f(a) - f(x));
                n++;
            }
            return x0;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static Complex secan(IComplex f, Complex a, Complex b, double eps = 1e-8)
        {
            Complex x1 = a;
            Complex x2 = b;
            Complex fb = f(b);
            Complex mpoint;
            int n = 0;

            while (Maths.Abs(f(x2)) > eps && n < short.MaxValue)
            {
                mpoint = x2 - (x2 - x1) * fb / (fb - f(x1));
                x1 = x2;
                x2 = mpoint;
                fb = f(x2);
                n++;
            }
            return x2;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static Complex falpo(IComplex f, Complex a, Complex b, double eps = 1e-8)
        {
            Complex x1 = a;
            Complex x2 = b;
            Complex fb = f(b);
            int n = 0;

            while (Maths.Abs(x2 - x1) > eps && n < short.MaxValue)
            {
                Complex xpoint = x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
                Complex fxpoint = f(xpoint);
                double s = fb.Real * fxpoint.Real;

                // sign
                if (s > 0)
                    x2 = xpoint;
                else
                    x1 = xpoint;

                if (Maths.Abs(fxpoint) < eps)
                    break;
                n++;
            }
            return x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
        }
        #endregion

        #region Enums
        /// <summary>
        /// Method for solving a nonlinear equation.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Bisection method.
            /// </summary>
            Bisection,
            /// <summary>
            /// Chord method.
            /// </summary>
            Chord,
            /// <summary>
            /// Secant method.
            /// </summary>
            Secant,
            /// <summary>
            /// False position method.
            /// </summary>
            FalsePosition,
            #endregion
        }
        #endregion
    }
}
