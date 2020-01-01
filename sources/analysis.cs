// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using UMapx.Core;
using UMapx.Decomposition;

namespace UMapx.Analysis
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                                 UMAPX.CORE
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Nonlinear solution
    /// <summary>
    /// Defines a class that implements the solution of a nonlinear equation.
    /// <remarks>
    /// This class is a solution to the problem of finding the root of a nonlinear equation of the form F(x) = 0.
    /// </remarks>
    /// </summary>
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
    #endregion

    #region Optimization methods
    /// <summary>
    /// Defines a class that implements an extremum search.
    /// <remarks>
    /// This class is a solution to the problem of finding the maximum and minimum points of the function F(x).
    /// </remarks>
    /// </summary>
    public class Optimization
    {
        #region Private data
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements an extremum search.
        /// </summary>
        /// <param name="eps">Epsilon [0, 1]</param>
        public Optimization(double eps = 1e-8)
        {
            this.Eps = eps;
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
        /// Returns the corresponding minimum of the function on the segment.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Start of line</param>
        /// <param name="b">End of line</param>
        /// <param name="max">Search maximum or minimum</param>
        /// <returns>Double precision floating point number</returns>
        public double Compute(IDouble function, double a, double b, bool max = false)
        {
            // max or min
            return (max) ? goldenMax(function, a, b, this.eps) : goldenMin(function, a, b, this.eps);
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
        private static double goldenMin(IDouble f, double a, double b, double eps = 1e-8)
        {
            double x1, x2;

            for (int i = 0; i < short.MaxValue; i++)
            {
                x1 = b - (b - a) / Maths.Phi;
                x2 = a + (b - a) / Maths.Phi;

                if (f(x1) > f(x2))
                    a = x1;
                else
                    b = x2;
                if (Math.Abs(b - a) < eps)
                    break;
            }
            return (a + b) / 2;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static double goldenMax(IDouble f, double a, double b, double eps = 1e-8)
        {
            double x1, x2;

            for (int i = 0; i < short.MaxValue; i++)
            {
                x1 = b - (b - a) / Maths.Phi;
                x2 = a + (b - a) / Maths.Phi;

                if (f(x1) < f(x2))
                    a = x1;
                else
                    b = x2;
                if (Math.Abs(b - a) < eps)
                    break;
            }
            return (a + b) / 2;
        }
        #endregion
    }
    #endregion

    #region Integral solution
    /// <summary>
    /// Defines a class that implements numerical integration.
    /// <remarks>
    /// This class is a solution to the problem of finding the value of the integral of the function F(x) within the values of a and b.
    /// </remarks>
    /// </summary>
    public class Integration
    {
        #region Private data
        private Integration.Method method;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements numerical integration.
        /// </summary>
        /// <param name="method">Integration method</param>
        public Integration(Integration.Method method = Method.Rectangle)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the integration method.
        /// </summary>
        public Integration.Method MethodType
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
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Double precision floating point number</returns>
        public double Compute(IDouble function, double a, double b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case Method.Midpoint:
                    return Integration.midp(function, a, b, n);

                case Method.Trapezoidal:
                    return Integration.trap(function, a, b, n);

                case Method.Simpson:
                    return Integration.simp(function, a, b, n);

                case Method.Romberg:
                    return Integration.romb(function, a, b, n);

                default:
                    return Integration.rect(function, a, b, n);
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="y">Function vector</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Double precision floating point number</returns>
        public double Compute(double[] y, double a, double b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case Method.Midpoint:
                    return Integration.midp(y, a, b, n);

                case Method.Trapezoidal:
                    return Integration.trap(y, a, b, n);

                case Method.Simpson:
                    return Integration.simp(y, a, b, n);

                default:
                    return Integration.rect(y, a, b, n);
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Complex number</returns>
        public Complex Compute(IComplex function, Complex a, Complex b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case Method.Midpoint:
                    return Integration.midp(function, a, b, n);

                case Method.Trapezoidal:
                    return Integration.trap(function, a, b, n);

                case Method.Simpson:
                    return Integration.simp(function, a, b, n);

                case Method.Romberg:
                    return Integration.romb(function, a, b, n);

                default:
                    return Integration.rect(function, a, b, n);
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="y">Function vector</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Complex number</returns>
        public Complex Compute(Complex[] y, Complex a, Complex b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case Method.Midpoint:
                    return Integration.midp(y, a, b, n);

                case Method.Trapezoidal:
                    return Integration.trap(y, a, b, n);

                case Method.Simpson:
                    return Integration.simp(y, a, b, n);

                default:
                    return Integration.rect(y, a, b, n);
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
        /// <param name="n"></param>
        /// <returns></returns>
        private static double rect(IDouble f, double a, double b, int n)
        {
            double sum = 0.0;
            double h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + i * h);
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static double rect(double[] y, double a, double b, int n)
        {
            double sum = 0.0;
            double h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * y[i];
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static double midp(IDouble f, double a, double b, int n)
        {
            // Midpoint
            double sum = 0.0;
            double h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + (i + 0.5) * h);
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static double midp(double[] y, double a, double b, int n)
        {
            double sum = 0.0;
            double h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += h * 0.5 * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static double trap(IDouble f, double a, double b, int n)
        {
            double sum = 0.0;
            double h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += 0.5 * h * (f(a + i * h) + f(a + (i + 1) * h));
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static double trap(double[] y, double a, double b, int n)
        {
            double sum = 0.0;
            double h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += 0.5 * h * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static double simp(IDouble f, double a, double b, int n)
        {
            if (n < 3) return double.NaN; //Need at least 3 points
            double sum = 0.0;
            double h = (b - a) / n;
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 1; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            else
            {
                sum = 3 * h * (f(a) + 3 * f(a + h) + 3 * f(a + 2 * h) + f(a + 3 * h)) / 8;
                for (int i = 3; i < n - 1; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static double simp(double[] y, double a, double b, int n)
        {
            double h = (b - a) / n;
            //Need at least 3 points
            if (n < 3 || h == 0) return double.NaN;
            double sum = 0.0;
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 1; i += 2)
                {
                    sum += h * (y[i] + 4 * y[i + 1] + y[i + 2]) / 3;
                }
            }
            else
            {
                sum = 3 * h * (y[0] + 3 * y[1] + 3 * y[2] + y[3]) / 8;
                for (int i = 3; i < n - 1; i += 2)
                {
                    sum += h * (y[i] + 4 * y[i + 1] + y[i + 2]) / 3;
                }
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="iterations"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static double romb(IDouble f, double a, double b, int iterations, double eps = 1e-8)
        {
            int n = 2;
            double h = b - a;
            double sum = 0.0;
            int j = 0;
            double[,] R = new double[iterations, iterations];
            R[1, 1] = h * (f(a) + f(b)) / 2.0;
            h = h / 2;
            R[2, 1] = R[1, 1] / 2 + h * f(a + h);
            R[2, 2] = (4 * R[2, 1] - R[1, 1]) / 3;
            for (j = 3; j <= iterations; j++)
            {
                n = 2 * n;
                h = h / 2;
                sum = 0.0;
                for (int k = 1; k <= n; k += 2)
                {
                    sum += f(a + k * h);
                }
                R[j, 1] = R[j - 1, 1] / 2 + h * sum;
                double factor = 4.0;
                for (int k = 2; k <= j; k++)
                {
                    R[j, k] = (factor * R[j, k - 1] - R[j - 1, k - 1]) / (factor - 1);
                    factor = factor * 4.0;
                }
                if (Math.Abs(R[j, j] - R[j, j - 1]) < eps * Math.Abs(R[j, j]))
                {
                    sum = R[j, j];
                    return sum;
                }
            }
            sum = R[n, n];
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex rect(IComplex f, Complex a, Complex b, int n)
        {
            Complex sum = 0.0;
            Complex h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + i * h);
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex rect(Complex[] y, Complex a, Complex b, int n)
        {
            Complex sum = 0.0;
            Complex h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * y[i];
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex midp(IComplex f, Complex a, Complex b, int n)
        {
            // Midpoint
            Complex sum = 0.0;
            Complex h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + (i + 0.5) * h);
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex midp(Complex[] y, Complex a, Complex b, int n)
        {
            Complex sum = 0.0;
            Complex h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += h * 0.5 * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex trap(IComplex f, Complex a, Complex b, int n)
        {
            Complex sum = 0.0;
            Complex h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += 0.5 * h * (f(a + i * h) + f(a + (i + 1) * h));
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex trap(Complex[] y, Complex a, Complex b, int n)
        {
            Complex sum = 0.0;
            Complex h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += 0.5 * h * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex simp(IComplex f, Complex a, Complex b, int n)
        {
            if (n < 3) return double.NaN; //Need at least 3 points
            Complex sum = 0.0;
            Complex h = (b - a) / n;
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 1; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            else
            {
                sum = 3 * h * (f(a) + 3 * f(a + h) + 3 * f(a + 2 * h) + f(a + 3 * h)) / 8;
                for (int i = 3; i < n - 1; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex simp(Complex[] y, Complex a, Complex b, int n)
        {
            Complex h = (b - a) / n;
            //Need at least 3 points
            if (n < 3 || h == 0) return double.NaN;
            Complex sum = 0.0;
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 1; i += 2)
                {
                    sum += h * (y[i] + 4 * y[i + 1] + y[i + 2]) / 3;
                }
            }
            else
            {
                sum = 3 * h * (y[0] + 3 * y[1] + 3 * y[2] + y[3]) / 8;
                for (int i = 3; i < n - 1; i += 2)
                {
                    sum += h * (y[i] + 4 * y[i + 1] + y[i + 2]) / 3;
                }
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="iterations"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static Complex romb(IComplex f, Complex a, Complex b, int iterations, double eps = 1e-8)
        {
            int n = 2;
            Complex h = b - a;
            Complex sum = 0.0;
            int j = 0;
            Complex[,] R = new Complex[iterations, iterations];
            R[1, 1] = h * (f(a) + f(b)) / 2.0;
            h = h / 2;
            R[2, 1] = R[1, 1] / 2 + h * f(a + h);
            R[2, 2] = (4 * R[2, 1] - R[1, 1]) / 3;
            for (j = 3; j <= iterations; j++)
            {
                n = 2 * n;
                h = h / 2;
                sum = 0.0;
                for (int k = 1; k <= n; k += 2)
                {
                    sum += f(a + k * h);
                }
                R[j, 1] = R[j - 1, 1] / 2 + h * sum;
                double factor = 4.0;
                for (int k = 2; k <= j; k++)
                {
                    R[j, k] = (factor * R[j, k - 1] - R[j - 1, k - 1]) / (factor - 1);
                    factor = factor * 4.0;
                }
                if (Maths.Abs(R[j, j] - R[j, j - 1]) < eps * Maths.Abs(R[j, j]))
                {
                    sum = R[j, j];
                    return sum;
                }
            }
            sum = R[n, n];
            return sum;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Integration method.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Rectangle method.
            /// </summary>
            Rectangle,
            /// <summary>
            /// Midpoint method.
            /// </summary>
            Midpoint,
            /// <summary>
            /// Trapezoidal method.
            /// </summary>
            Trapezoidal,
            /// <summary>
            /// Simpson method.
            /// </summary>
            Simpson,
            /// <summary>
            /// Romberg method.
            /// </summary>
            Romberg,
            #endregion
        }
        #endregion
    }
    #endregion

    #region Numeric differentiation
    /// <summary>
    /// Defines a class that implements numerical differentiation.
    /// </summary>
    public class Differentation
    {
        #region Private data
        private int points;
        #endregion

        #region Components
        /// <summary>
        /// Initializes a class that implements numerical differentiation.
        /// </summary>
        /// <param name="points">Number of interpolation points</param>
        public Differentation(int points)
        {
            this.Points = points;
        }
        /// <summary>
        /// Gets or sets the number of interpolation points.
        /// </summary>
        public int Points
        {
            get
            {
                return this.points;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.points = value;
            }
        }
        /// <summary>
        /// Returns the value of a derived function.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="x">Argument value</param>
        /// <param name="h">Step</param>
        /// <param name="order">Order</param>
        /// <returns>Double precision floating point number</returns>
        public double Compute(IDouble function, double x, double h, int order)
        {
            // exception
            if (order > this.points)
                throw new Exception("The order of the derivative cannot be greater than the number of interpolation points");
            if (order < 0)
                throw new Exception("The derivative order cannot be less than 0");

            // Create the interpolation points
            int length = this.points + 1;
            double[,] coefficients = Differentation.GetCoefficients(length);
            double sum = 0.0;

            // do job
            for (int i = 0, center = 0; i < length; i++)
            {
                sum += coefficients[order, i] * function(x + center * h);
                center++;
            }

            // result
            return sum / Math.Pow(h, order);
        }
        /// <summary>
        /// Returns the value of a derived function.
        /// </summary>
        /// <param name="y">Function vector</param>
        /// <param name="index">Index of argument</param>
        /// <param name="h">Step</param>
        /// <param name="order">Order</param>
        /// <returns>Double precision floating point number</returns>
        public double Compute(double[] y, int index, double h, int order)
        {
            // exception
            if (order > this.points)
                throw new Exception("The order of the derivative cannot be greater than the number of interpolation points");
            if (order < 0)
                throw new Exception("The derivative order cannot be less than 0");

            // Create the interpolation points
            int length = this.points + 1;
            double[,] coefficients = Differentation.GetCoefficients(length);
            double sum = 0.0;

            // do job
            for (int i = 0, center = 0; i < length; i++)
            {
                sum += coefficients[order, i] * y[index + i];
                center++;
            }

            // result
            return sum / Math.Pow(h, order);
        }
        /// <summary>
        /// Returns the value of a derived function.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="x">Argument value</param>
        /// <param name="h">Step</param>
        /// <param name="order">Order</param>
        /// <returns>Complex number</returns>
        public Complex Compute(IComplex function, Complex x, Complex h, int order)
        {
            // exception
            if (order > this.points)
                throw new Exception("The order of the derivative cannot be greater than the number of interpolation points");
            if (order < 0)
                throw new Exception("The derivative order cannot be less than 0");

            // Create the interpolation points
            int length = this.points + 1;
            double[,] coefficients = Differentation.GetCoefficients(length);
            Complex sum = 0.0;

            // do job
            for (int i = 0, center = 0; i < length; i++)
            {
                sum += coefficients[order, i] * function(x + center * h);
                center++;
            }

            // result
            return sum / Maths.Pow(h, order);
        }
        /// <summary>
        /// Returns the value of a derived function.
        /// </summary>
        /// <param name="y">Function vector</param>
        /// <param name="index">Index of argument</param>
        /// <param name="h">Step</param>
        /// <param name="order">Order</param>
        /// <returns>Complex number</returns>
        public Complex Compute(Complex[] y, int index, double h, int order)
        {
            // exception
            if (order > this.points)
                throw new Exception("The order of the derivative cannot be greater than the number of interpolation points");
            if (order < 0)
                throw new Exception("The derivative order cannot be less than 0");

            // Create the interpolation points
            int length = this.points + 1;
            double[,] coefficients = Differentation.GetCoefficients(length);
            Complex sum = 0.0;

            // do job
            for (int i = 0, center = 0; i < length; i++)
            {
                sum += coefficients[order, i] * y[index + i];
                center++;
            }

            // result
            return sum / Math.Pow(h, order);
        }
        #endregion

        #region Static voids
        /// <summary>
        /// Returns the matrix of interpolation coefficients.
        /// </summary>
        /// <param name="points">Number of points</param>
        /// <returns>Matrix</returns>
        public static double[,] GetCoefficients(int points)
        {
            // Compute difference coefficient table
            double fac = Special.Factorial(points);
            double[,] deltas = new double[points, points];
            double h, delta;
            int j, k;

            // do job
            for (j = 0; j < points; j++)
            {
                h = 1.0;
                delta = j;

                for (k = 0; k < points; k++)
                {
                    deltas[j, k] = h / Special.Factorial(k);
                    h *= delta;
                }
            }

            // matrix invert
            deltas = Matrice.Invert(deltas);

            //// rounding
            //for (j = 0; j < points; j++)
            //    for (k = 0; k < points; k++)
            //        deltas[j, k] = (Math.Round(deltas[j, k] * fac, MidpointRounding.AwayFromZero)) / fac;

            return deltas;
        }
        #endregion
    }
    #endregion

    #region Differential equation solution
    /// <summary>
    /// Defines a class that implements a solution to a differential equation.
    /// <remarks>
    /// This class is a solution to the Cauchy problem for the ordinary differential equation y' = F(x, y).
    /// </remarks>
    /// </summary>
    public class Differential
    {
        #region Private data
        private Differential.Method method;
        #endregion

        #region Diferentiation components
        /// <summary>
        /// Initializes a class that implements the solution of a differential equation.
        /// </summary>
        /// <param name="method">Differentiation method</param>
        public Differential(Differential.Method method = Method.RungeKutta4)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the differentiation method.
        /// </summary>
        public Differential.Method MethodType
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
        /// Returns the value of a differential equation.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <returns>Array of function values</returns>
        public double[] Compute(IDoubleMesh function, double[] x, double y0)
        {
            // chose method of differentiation
            switch (method)
            {
                case Method.Euler:
                    return Differential.euler(function, x, y0);

                case Method.Fehlberg:
                    return Differential.fehlberg(function, x, y0);

                case Method.RungeKutta4:
                    return Differential.rungeKutta4(function, x, y0);

                default:
                    return Differential.rungeKutta2(function, x, y0);
            }
        }
        /// <summary>
        /// Returns the value of a differential equation.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <returns>Array of function values</returns>
        public Complex[] Compute(IComplexMesh function, Complex[] x, Complex y0)
        {
            // chose method of differentiation
            switch (method)
            {
                case Method.Euler:
                    return Differential.euler(function, x, y0);

                case Method.Fehlberg:
                    return Differential.fehlberg(function, x, y0);

                case Method.RungeKutta4:
                    return Differential.rungeKutta4(function, x, y0);

                default:
                    return Differential.rungeKutta2(function, x, y0);
            }
        }
        #endregion

        #region Recompute voids
        /// <summary>
        /// Returns the value of a differential equation calculated by the Adams-Bashfort method.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <param name="order">Order</param>
        /// <returns>Array of function values</returns>
        public double[] Compute(IDoubleMesh function, double[] x, double y0, int order = 2)
        {
            int n = x.Length - 1;

            // if order more than 1
            // Adams-Bashfort method
            if (order > 1 && order < n)
            {
                // params
                int i, j, k = order + 1;
                double[] y = new double[n];
                double[] r = new double[k];
                double[] c = Differential.GetCoefficients(order);
                double h, t, sum;

                // compute first points by order
                for (i = 0; i < k; i++)
                    r[i] = x[i];

                // classic differential
                r = this.Compute(function, r, y0);

                for (i = 0; i < order; i++)
                    y[i] = r[i];

                // Adams-Bashforth method
                // for order
                for (i = order; i < n; i++)
                {
                    sum = y[i - 1];

                    for (j = 0; j < order; j++)
                    {
                        t = x[i - j];
                        h = t - x[i - j - 1];
                        sum += h * c[j] * function(t, y[i - j - 1]);
                    }

                    y[i] = sum;
                }

                return y;
            }

            // classic differential
            return this.Compute(function, x, y0);
        }
        /// <summary>
        /// Returns the value of a differential equation calculated by the Adams-Bashfort method.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <param name="order">Order</param>
        /// <returns>Array of function values</returns>
        public Complex[] Compute(IComplexMesh function, Complex[] x, Complex y0, int order = 2)
        {
            int n = x.Length - 1;

            // if order more than 1
            // Adams-Bashfort method
            if (order > 1 && order < n)
            {
                // params
                int i, j, k = order + 1;
                Complex[] y = new Complex[n];
                Complex[] r = new Complex[k];
                double[] c = Differential.GetCoefficients(order);
                Complex h, t, sum;

                // compute first points by order
                for (i = 0; i < k; i++)
                    r[i] = x[i];

                // classic differential
                r = this.Compute(function, r, y0);

                for (i = 0; i < order; i++)
                    y[i] = r[i];

                // Adams-Bashforth method
                // for order
                for (i = order; i < n; i++)
                {
                    sum = y[i - 1];

                    for (j = 0; j < order; j++)
                    {
                        t = x[i - j];
                        h = t - x[i - j - 1];
                        sum += h * c[j] * function(t, y[i - j - 1]);
                    }

                    y[i] = sum;
                }

                return y;
            }

            // classic differential
            return this.Compute(function, x, y0);
        }
        #endregion

        #region Adams-Bashforth
        /// <summary>
        /// Returns an array of coefficient values for the Adams-Bashfort formula.
        /// </summary>
        /// <param name="order">Order</param>
        /// <returns>Array</returns>
        public static double[] GetCoefficients(int order)
        {
            double[,] A = new double[order, order];
            double[] c = new double[order];
            int i, j;

            for (i = 0; i < order; i++)
            {
                for (j = 0; j < order; j++)
                {
                    A[i, j] = Math.Pow(j, i);
                }
                c[i] = Math.Pow(-1, i) / (i + 1);
            }

            return A.Solve(c);
            //return c.Dot(A.Invert());
        }
        #endregion

        #region Runge-Kutta
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static double[] euler(IDoubleMesh f, double[] x, double y0)
        {
            int n = x.Length - 1;
            double xnew, ynew = y0, h;
            double[] result = new double[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                ynew = ynew + f(xnew, ynew) * h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static double[] rungeKutta2(IDoubleMesh f, double[] x, double y0)
        {
            int n = x.Length - 1;
            double xnew, ynew = y0, h, k1, k2;
            double[] result = new double[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.5 * h, ynew + 0.5 * k1);
                ynew = ynew + k2;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static double[] rungeKutta4(IDoubleMesh f, double[] x, double y0)
        {
            int n = x.Length - 1;
            double xnew, ynew = y0, h, k1, k2, k3, k4;
            double[] result = new double[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.5 * h, ynew + 0.5 * k1);
                k3 = h * f(xnew + 0.5 * h, ynew + 0.5 * k2);
                k4 = h * f(xnew + h, ynew + k3);
                ynew = ynew + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static double[] fehlberg(IDoubleMesh f, double[] x, double y0)
        {
            int n = x.Length - 1;
            double xnew, ynew = y0, h, k1, k2, k3, k4, k5, k6;
            double[] result = new double[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.25 * h, ynew + 0.25 * k1);
                k3 = h * f(xnew + 3 * h / 8, ynew + 3 * k1 / 32 + 9 * k2 / 32);
                k4 = h * f(xnew + 12 * h / 13, ynew + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197);
                k5 = h * f(xnew + h, ynew + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104);
                k6 = h * f(xnew + 0.5 * h, ynew - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40);
                ynew = ynew + 25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - 0.2 * k5;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static Complex[] euler(IComplexMesh f, Complex[] x, Complex y0)
        {
            int n = x.Length - 1;
            Complex xnew, ynew = y0, h;
            Complex[] result = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                ynew = ynew + f(xnew, ynew) * h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static Complex[] rungeKutta2(IComplexMesh f, Complex[] x, Complex y0)
        {
            int n = x.Length - 1;
            Complex xnew, ynew = y0, h, k1, k2;
            Complex[] result = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.5 * h, ynew + 0.5 * k1);
                ynew = ynew + k2;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static Complex[] rungeKutta4(IComplexMesh f, Complex[] x, Complex y0)
        {
            int n = x.Length - 1;
            Complex xnew, ynew = y0, h, k1, k2, k3, k4;
            Complex[] result = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.5 * h, ynew + 0.5 * k1);
                k3 = h * f(xnew + 0.5 * h, ynew + 0.5 * k2);
                k4 = h * f(xnew + h, ynew + k3);
                ynew = ynew + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static Complex[] fehlberg(IComplexMesh f, Complex[] x, Complex y0)
        {
            int n = x.Length - 1;
            Complex xnew, ynew = y0, h, k1, k2, k3, k4, k5, k6;
            Complex[] result = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.25 * h, ynew + 0.25 * k1);
                k3 = h * f(xnew + 3 * h / 8, ynew + 3 * k1 / 32 + 9 * k2 / 32);
                k4 = h * f(xnew + 12 * h / 13, ynew + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197);
                k5 = h * f(xnew + h, ynew + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104);
                k6 = h * f(xnew + 0.5 * h, ynew - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40);
                ynew = ynew + 25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - 0.2 * k5;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Differentiation method
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Euler method.
            /// </summary>
            Euler,
            /// <summary>
            /// The second-order Runge-Kutta method.
            /// </summary>
            RungeKutta2,
            /// <summary>
            /// Fourth-order Runge-Kutta method.
            /// </summary>
            RungeKutta4,
            /// <summary>
            /// Felberg's method.
            /// </summary>
            Fehlberg,
            #endregion
        }
        #endregion
    }
    #endregion

    #region Interpolation methods
    /// <summary>
    /// Defines a class that implements interpolation.
    /// <remarks>
    /// This class is a solution to the problem of finding an intermediate value of the function F(x).
    /// </remarks>
    /// </summary>
    public class Interpolation
    {
        #region Private data
        private Interpolation.Method method;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements interpolation.
        /// </summary>
        /// <param name="method">Interpolation method</param>
        public Interpolation(Interpolation.Method method = Method.Lagrange)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the interpolation method.
        /// </summary>
        public Interpolation.Method MethodType
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
        /// Returns the value of a function at a point.
        /// <remarks>
        /// In this case, only bilinear interpolation is used.
        /// </remarks>
        /// </summary>
        /// <param name="x">Array of values of the first argument</param>
        /// <param name="y">Array of values of the second argument</param>
        /// <param name="z">Function matrix</param>
        /// <param name="xl">The value of the first argument to calculate</param>
        /// <param name="yl">The value of the second argument to calculate</param>
        /// <returns>Double precision floating point number</returns>
        public double Compute(double[] x, double[] y, double[,] z, double xl, double yl)
        {
            return bilinear(x, y, z, xl, yl);
        }
        /// <summary>
        /// Returns the value of a function at a point.
        /// </summary>
        /// <param name="x">Array of values of the argument</param>
        /// <param name="y">Array of values of the function</param>
        /// <param name="xl">The value of the argument to calculate</param>
        /// <returns>Double precision floating point number</returns>
        public double Compute(double[] x, double[] y, double xl)
        {
            // chose method of interpolation
            switch (method)
            {
                case Method.Lagrange:
                    return Interpolation.lagra(x, y, xl);

                case Method.Newton:
                    return Interpolation.newto(x, y, xl);

                case Method.Barycentric:
                    return Interpolation.baryc(x, y, xl);

                default:
                    return Interpolation.linear(x, y, xl);
            }
        }
        /// <summary>
        /// Returns the value of a function at a point.
        /// </summary>
        /// <param name="x">Array of values of the argument</param>
        /// <param name="y">Array of values of the function</param>
        /// <param name="xl">The value of the argument to calculate</param>
        /// <returns>Complex number</returns>
        public Complex Compute(Complex[] x, Complex[] y, Complex xl)
        {
            // chose method of interpolation
            switch (method)
            {
                case Method.Newton:
                    return Interpolation.newto(x, y, xl);

                case Method.Barycentric:
                    return Interpolation.baryc(x, y, xl);

                default:
                    return Interpolation.lagra(x, y, xl);
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xl"></param>
        /// <returns></returns>
        private static double linear(double[] x, double[] y, double xl)
        {
            double yval = 0.0;
            int length = x.Length - 1;

            for (int i = 0; i < length; i++)
            {
                if (xl >= x[i] && xl < x[i + 1])
                {
                    yval = y[i] + (xl - x[i]) * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
                }
            }
            return yval;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <returns></returns>
        private static double bilinear(double[] x, double[] y, double[,] z, double xval, double yval)
        {
            double zval = 0.0;
            int xlength = x.Length - 1;
            int ylength = y.Length - 1;

            for (int i = 0; i < xlength; i++)
            {
                for (int j = 0; j < ylength; j++)
                {
                    if (xval >= x[i] && xval < x[i + 1] && yval >= y[j] && yval < y[j + 1])
                    {
                        zval = z[i, j] * (x[i + 1] - xval) * (y[j + 1] - yval) / (x[i + 1] - x[i]) / (y[j + 1] - y[j]) +
                        z[i + 1, j] * (xval - x[i]) * (y[j + 1] - yval) / (x[i + 1] - x[i]) / (y[j + 1] - y[j]) +
                        z[i, j + 1] * (x[i + 1] - xval) * (yval - y[j]) / (x[i + 1] - x[i]) / (y[j + 1] - y[j]) +
                        z[i + 1, j + 1] * (xval - x[i]) * (yval - y[j]) / (x[i + 1] - x[i]) / (y[j + 1] - y[j]);
                    }
                }
            }
            return zval;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static double lagra(double[] x, double[] y, double xval)
        {
            double yval = 0.0;
            double Products = y[0];
            int length = x.Length;
            int i, j;

            for (i = 0; i < length; i++)
            {
                Products = y[i];
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        Products *= (xval - x[j]) / (x[i] - x[j]);
                    }
                }
                yval += Products;
            }
            return yval;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static double newto(double[] x, double[] y, double xval)
        {
            double yval;
            int length = x.Length;
            double[] tarray = new double[length];
            int i, j;

            for (i = 0; i < length; i++)
            {
                tarray[i] = y[i];
            }
            for (i = 0; i < length - 1; i++)
            {
                for (j = length - 1; j > i; j--)
                {
                    tarray[j] = (tarray[j - 1] - tarray[j]) / (x[j - 1 - i] - x[j]);
                }
            }
            yval = tarray[length - 1];
            for (i = length - 2; i >= 0; i--)
            {
                yval = tarray[i] + (xval - x[i]) * yval;
            }
            return yval;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static double baryc(double[] x, double[] y, double xval)
        {
            double product;
            double deltaX;
            double bc1 = 0;
            double bc2 = 0;
            int length = x.Length;
            double[] weights = new double[length];
            int i, j;

            for (i = 0; i < length; i++)
            {
                product = 1;
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        product *= (x[i] - x[j]);
                        weights[i] = 1.0 / product;
                    }
                }
            }

            for (i = 0; i < length; i++)
            {
                deltaX = weights[i] / (xval - x[i]);
                bc1 += y[i] * deltaX;
                bc2 += deltaX;
            }
            return bc1 / bc2;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static Complex lagra(Complex[] x, Complex[] y, Complex xval)
        {
            Complex yval = 0.0;
            Complex Products = y[0];
            int length = x.Length;
            int i, j;

            for (i = 0; i < length; i++)
            {
                Products = y[i];
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        Products *= (xval - x[j]) / (x[i] - x[j]);
                    }
                }
                yval += Products;
            }
            return yval;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static Complex newto(Complex[] x, Complex[] y, Complex xval)
        {
            Complex yval;
            int length = x.Length;
            Complex[] tarray = new Complex[length];
            int i, j;

            for (i = 0; i < length; i++)
            {
                tarray[i] = y[i];
            }
            for (i = 0; i < length - 1; i++)
            {
                for (j = length - 1; j > i; j--)
                {
                    tarray[j] = (tarray[j - 1] - tarray[j]) / (x[j - 1 - i] - x[j]);
                }
            }
            yval = tarray[length - 1];
            for (i = length - 2; i >= 0; i--)
            {
                yval = tarray[i] + (xval - x[i]) * yval;
            }
            return yval;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static Complex baryc(Complex[] x, Complex[] y, Complex xval)
        {
            Complex product;
            Complex deltaX;
            Complex bc1 = 0;
            Complex bc2 = 0;
            int length = x.Length;
            Complex[] weights = new Complex[length];
            int i, j;

            for (i = 0; i < length; i++)
            {
                product = 1;
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        product *= (x[i] - x[j]);
                        weights[i] = 1.0 / product;
                    }
                }
            }

            for (i = 0; i < length; i++)
            {
                deltaX = weights[i] / (xval - x[i]);
                bc1 += y[i] * deltaX;
                bc2 += deltaX;
            }
            return bc1 / bc2;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Interpolation method.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Linear method.
            /// </summary>
            Linear,
            /// <summary>
            /// Lagrange's method.
            /// </summary>
            Lagrange,
            /// <summary>
            /// Newton's method.
            /// </summary>
            Newton,
            /// <summary>
            /// Barycentric method.
            /// </summary>
            Barycentric,
            #endregion
        }
        #endregion
    }
    #endregion

    #region Approximation methods
    /// <summary>
    /// Defines the least squares approximation class.
    /// <remarks>
    /// This class is a solution to the problem of finding the function A (x) ≈ F (x), where F (x) is the original function.
    /// More information can be found on the website:
    /// http://simenergy.ru/math-analysis/digital-processing/85-ordinary_least_squares
    /// </remarks>
    /// </summary>
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
            double[] cf = null;
            double error = 0;
            string equation = null;

            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, ref cf, ref error, ref equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, ref cf, ref error, ref equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, ref cf, ref error, ref equation);

                default:
                    return Approximation.powr(x, y, power, ref cf, ref error, ref equation);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <returns>Array</returns>
        public double[] Compute(double[] x, double[] y, ref double[] cf)
        {
            double error = 0;
            string equation = null;

            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, ref cf, ref error, ref equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, ref cf, ref error, ref equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, ref cf, ref error, ref equation);

                default:
                    return Approximation.powr(x, y, power, ref cf, ref error, ref equation);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="error">Error</param>
        /// <returns>Array</returns>
        public double[] Compute(double[] x, double[] y, ref double[] cf, ref double error)
        {
            string equation = null;

            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, ref cf, ref error, ref equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, ref cf, ref error, ref equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, ref cf, ref error, ref equation);

                default:
                    return Approximation.powr(x, y, power, ref cf, ref error, ref equation);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="error">Error</param>
        /// <param name="equation">Equation</param>
        /// <returns>Array</returns>
        public double[] Compute(double[] x, double[] y, ref double[] cf, ref double error, ref string equation)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, ref cf, ref error, ref equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, ref cf, ref error, ref equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, ref cf, ref error, ref equation);

                default:
                    return Approximation.powr(x, y, power, ref cf, ref error, ref equation);
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
            Complex[] cf = null;
            Complex error = 0;
            string equation = null;

            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, ref cf, ref error, ref equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, ref cf, ref error, ref equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, ref cf, ref error, ref equation);

                default:
                    return Approximation.powr(x, y, power, ref cf, ref error, ref equation);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <returns>Array</returns>
        public Complex[] Compute(Complex[] x, Complex[] y, ref Complex[] cf)
        {
            Complex error = 0;
            string equation = null;

            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, ref cf, ref error, ref equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, ref cf, ref error, ref equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, ref cf, ref error, ref equation);

                default:
                    return Approximation.powr(x, y, power, ref cf, ref error, ref equation);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="error">Error</param>
        /// <returns>Array</returns>
        public Complex[] Compute(Complex[] x, Complex[] y, ref Complex[] cf, ref Complex error)
        {
            string equation = null;

            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, ref cf, ref error, ref equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, ref cf, ref error, ref equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, ref cf, ref error, ref equation);

                default:
                    return Approximation.powr(x, y, power, ref cf, ref error, ref equation);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="error">Error</param>
        /// <param name="equation">Equation</param>
        /// <returns>Array</returns>
        public Complex[] Compute(Complex[] x, Complex[] y, ref Complex[] cf, ref Complex error, ref string equation)
        {
            // chose method of approximation
            switch (method)
            {
                case Method.Polynomial:
                    return Approximation.poly(x, y, power, ref cf, ref error, ref equation);

                case Method.Logarithmic:
                    return Approximation.logc(x, y, power, ref cf, ref error, ref equation);

                case Method.Exponential:
                    return Approximation.expn(x, y, power, ref cf, ref error, ref equation);

                default:
                    return Approximation.powr(x, y, power, ref cf, ref error, ref equation);
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
        private static double[] poly(double[] x, double[] y, int power, ref double[] cf, ref double error, ref string equation)
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
        private static Complex[] poly(Complex[] x, Complex[] y, int power, ref Complex[] cf, ref Complex error, ref string equation)
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
        private static double[] logc(double[] x, double[] y, int power, ref double[] cf, ref double error, ref string equation)
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
        private static Complex[] logc(Complex[] x, Complex[] y, int power, ref Complex[] cf, ref Complex error, ref string equation)
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
        private static double[] expn(double[] x, double[] y, int power, ref double[] cf, ref double error, ref string equation)
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
        private static Complex[] expn(Complex[] x, Complex[] y, int power, ref Complex[] cf, ref Complex error, ref string equation)
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
        private static double[] powr(double[] x, double[] y, int power, ref double[] cf, ref double error, ref string equation)
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
        private static Complex[] powr(Complex[] x, Complex[] y, int power, ref Complex[] cf, ref Complex error, ref string equation)
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

        #region Complex components
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
    #endregion

    #region Roots solution
    /// <summary>
    /// Defines a class for solving equations using the spectral decomposition of a matrix.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.mathworks.com/help/matlab/ref/roots.html
    /// </remarks>
    /// </summary>
    public class Roots
    {
        #region Private data
        private EVD eig;
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class of equations using the spectral decomposition of a matrix.
        /// </summary>
        /// <param name="eps">Epsilon [0, 1]</param>
        public Roots(double eps = 1e-16)
        {
            this.Eps = eps;
        }
        /// <summary>
        /// Gets or sets an error [0, 1].
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
        /// Returns a column vector corresponding to the numerical solution of the polynomial: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="polynomial">Polynomial</param>
        /// <returns>Array</returns>
        public Complex[] Compute(double[] polynomial)
        {
            // MATLAB roots method
            // represented by Asiryan Valeriy, 2018.
            // properties of polynomial:
            int length = polynomial.Length;
            int i, index = -1;

            // finding non-zero element:
            for (i = 0; i < length; i++)
            {
                if (polynomial[i] != 0)
                {
                    index = i;
                    break;
                }
            }

            // return null array:
            if (index == -1)
            {
                return new Complex[0];
            }

            // get scaling factor:
            int m = length - index - 1;
            double scale = polynomial[index];
            double[] c = new double[m];

            // create new polynomial:
            for (i = 0; i < m; i++)
            {
                c[i] = polynomial[i + index + 1] / scale;
            }
            
            // Eigen-value decomposition for
            // companion matrix:
            eig = new EVD(Matrice.Companion(c), this.eps);

            // Complex result:
            return eig.D;
        }
        /// <summary>
        /// Returns a column vector of polynomial coefficients: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="roots">Roots</param>
        /// <returns>Array</returns>
        public double[] Compute(Complex[] roots)
        {
            // MATLAB roots method
            // represented by Asiryan Valeriy, 2018.
            // properties of polynomial:
            int length = roots.Length, m = length + 1, j, i;

            // arrays:
            Complex[] v = new Complex[length];
            Complex[] p = new Complex[m];

            // point:
            p[0] = 1.0;

            // create new polynomial:
            for (j = 0; j < length; j++)
            {
                // right part:
                for (i = 0; i <= j; i++)
                {
                    v[i] = roots[j] * p[i];
                }
                // left part:
                for (i = 0; i <= j; i++)
                {
                    p[i + 1] -= v[i];
                }
            }

            // Real result:
            return p.Real();
        }
        #endregion
    }
    #endregion
}
