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

namespace UMapx.Core
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
    /// Определяет класс, реализующий решение нелинейного уравнения.
    /// </summary>
    public class Nonlinear
    {
        #region Private data
        private Nonlinear.Method method;
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий решение нелинейного уравнения.
        /// </summary>
        /// <param name="eps">Погрещность [0, 1]</param>
        /// <param name="method">Метод решения нелинейного уравнения</param>
        public Nonlinear(double eps = 1e-8, Nonlinear.Method method = Method.Bisection)
        {
            this.method = method;
            this.Eps = eps;
        }
        /// <summary>
        /// Получает или задает метод решения нелинейного уравнения.
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
        /// Получает или задает значение погрешности [0, 1].
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
        /// Получает значения корня нелинейного уравнения.
        /// </summary>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <param name="a">Начало отрезка</param>
        /// <param name="b">Конец отрезка</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(IDouble function, double a, double b)
        {
            // chose method of nonlinear
            switch (method)
            {
                case Method.Chord:
                    return Nonlinear.chord(function, a, b, this.eps);
                case Method.FalsePosition:
                    return Nonlinear.chord(function, a, b, this.eps);
                case Method.Secant:
                    return Nonlinear.chord(function, a, b, this.eps);

                default:
                    return Nonlinear.bisec(function, a, b, this.eps);
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

            while (Maths.Abs(f(x0) / b) > eps && n < short.MaxValue)
            {
                x = x0;
                x0 = x - (f(x) * (a - x)) / (f(a) - f(x));
                n++;
            }
            return x0;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Метод решения нелинейного уравнения.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Метод половинного деления.
            /// </summary>
            Bisection,
            /// <summary>
            /// Метод Ньютона.
            /// </summary>
            Chord,
            /// <summary>
            /// Метод секансов.
            /// </summary>
            Secant,
            /// <summary>
            /// Метод ложной точки.
            /// </summary>
            FalsePosition,
            #endregion
        }
        #endregion
    }
    #endregion

    #region Interpolation methods
    /// <summary>
    /// Определяет класс, реализующий интерполяцию.
    /// </summary>
    public class Interpolation
    {
        #region Private data
        private Interpolation.Method method;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий интерполяцию.
        /// </summary>
        /// <param name="method">Метод интерполяции</param>
        public Interpolation(Interpolation.Method method = Method.Lagrange)
        {
            this.method = method;
        }
        /// <summary>
        /// Получает или задает метод интерполяции.
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
        /// Возвращает значение функции в точке.
        /// <remarks>
        /// В данном случае используется только биленейная интерполяция.
        /// </remarks>
        /// </summary>
        /// <param name="x">Массив значений первого аргумента</param>
        /// <param name="y">Массив значений второго аргумента</param>
        /// <param name="z">Матрица значений функции</param>
        /// <param name="xl">Значение первого аргумента для вычисления</param>
        /// <param name="yl">Значение второго аргумента для вычисления</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(double[] x, double[] y, double[,] z, double xl, double yl)
        {
            return bilinear(x, y, z, xl, yl);
        }
        /// <summary>
        /// Возвращает значение функции в точке.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="xl">Значение аргумента для вычисления</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        public static double lagra(double[] x, double[] y, double xval)
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
        public static double newto(double[] x, double[] y, double xval)
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
        public static double baryc(double[] x, double[] y, double xval)
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
        #endregion

        #region Enums
        /// <summary>
        /// Метод интерполяции.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Линейный метод.
            /// </summary>
            Linear,
            /// <summary>
            /// Метод Лагранжа.
            /// </summary>
            Lagrange,
            /// <summary>
            /// Метод Ньютона.
            /// </summary>
            Newton,
            /// <summary>
            /// Барицентрический метод.
            /// </summary>
            Barycentric,
            #endregion
        }
        #endregion
    }
    #endregion

    #region Optimization methods
    /// <summary>
    /// Определяет класс, реализующий поиск экстремума.
    /// </summary>
    public class Optimization
    {
        #region Private data
        /// <summary>
        /// Погрешность вычислений.
        /// </summary>
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий поиск экстремума.
        /// </summary>
        /// <param name="eps">Погрешность [0, 1]</param>
        public Optimization(double eps = 1e-8)
        {
            this.Eps = eps;
        }
        /// <summary>
        /// Получает или задает значение погрешности [0, 1].
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
        /// Возвращает соответствующий минимум функции на отреке.
        /// </summary>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <param name="a">Начало отрезка</param>
        /// <param name="b">Конец отрезка</param>
        /// <param name="max">Искать максимум или минимум</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
    /// Определяет класс, реализующий численное интегрирование.
    /// </summary>
    public class Integration
    {
        #region Private data
        private Integration.Method method;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий численное интегрирование.
        /// </summary>
        /// <param name="method">Метод интегрирования</param>
        public Integration(Integration.Method method = Method.Rectangle)
        {
            this.method = method;
        }
        /// <summary>
        /// Получает или задает метод интегрирования.
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
        /// Возвращает значение интеграла функции.
        /// </summary>
        /// <param name="function">Делегат непрерывной фукнции</param>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <param name="n">Количество разбиений</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает значение интеграла функции.
        /// </summary>
        /// <param name="y">Вектор значений фукнции</param>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <param name="n">Количество разбиений</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        #endregion

        #region Enums
        /// <summary>
        /// Метод интегрирования.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Метод прямоугольников.
            /// </summary>
            Rectangle,
            /// <summary>
            /// Метод средней точки.
            /// </summary>
            Midpoint,
            /// <summary>
            /// Метод трапеций.
            /// </summary>
            Trapezoidal,
            /// <summary>
            /// Метод Симпсона.
            /// </summary>
            Simpson,
            /// <summary>
            /// Метод Ромберга.
            /// </summary>
            Romberg,
            #endregion
        }
        #endregion
    }
    #endregion

    #region Differential solution
    /// <summary>
    /// Определяет класс, реализующий решение дифференциального уравнения.
    /// </summary>
    public class Diferentiation
    {
        #region Private data
        private Diferentiation.Method method;
        #endregion

        #region Diferentiation components
        /// <summary>
        /// Инициализирует класс, реализующий решение дифференциального уравнения.
        /// </summary>
        /// <param name="method">Метод дифференцирования</param>
        public Diferentiation(Diferentiation.Method method = Method.RungeKutta4)
        {
            this.method = method;
        }
        /// <summary>
        /// Получает или задает метод дифференцирования.
        /// </summary>
        public Diferentiation.Method MethodType
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
        /// Возвращает значение дифференциального уравнения.
        /// </summary>
        /// <param name="function">Делегат непрерывной функции, зависящей от двух переменных</param>
        /// <param name="x0">Начало отрезка</param>
        /// <param name="x1">Конец отрезка</param>
        /// <param name="h">Шаг</param>
        /// <param name="y0">Значение</param>
        /// <returns>Значение функции</returns>
        public double Compute(IDoubleMesh function, double x0, double x1, double h, double y0)
        {
            // chose method of differentiation
            switch (method)
            {
                case Method.Euler:
                    return Diferentiation.euler(function, x0, x1, h, y0);

                case Method.Fehlberg:
                    return Diferentiation.rungeKuttaFehlberg(function, x0, x1, h, y0, 1e-8);

                case Method.RungeKutta4:
                    return Diferentiation.rungeKutta4(function, x0, x1, h, y0);

                default:
                    return Diferentiation.rungeKutta2(function, x0, x1, h, y0);
            }
        }
        #endregion

        #region Private double voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x0"></param>
        /// <param name="x"></param>
        /// <param name="h"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static double euler(IDoubleMesh f, double x0, double x, double h, double y0)
        {
            int n = 0;
            double xnew, ynew, result = double.NaN;
            if (x <= x0)
                result = y0;
            else if (x > x0)
            {
                do
                {
                    if (h > x - x0) h = x - x0;
                    ynew = y0 + f(x0, y0) * h;
                    xnew = x0 + h;
                    x0 = xnew;
                    y0 = ynew;
                    n++;
                } while (x0 < x && n < short.MaxValue);
                result = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x0"></param>
        /// <param name="x"></param>
        /// <param name="h"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static double rungeKutta2(IDoubleMesh f, double x0, double x, double h, double y0)
        {
            int n = 0;
            double xnew, ynew, k1, k2, result = double.NaN;
            if (x == x0)
                result = y0;
            else if (x > x0)
            {
                do
                {
                    if (h > x - x0) h = x - x0;
                    k1 = h * f(x0, y0);
                    k2 = h * f(x0 + 0.5 * h, y0 + 0.5 * k1);
                    ynew = y0 + k2;
                    xnew = x0 + h;
                    x0 = xnew;
                    y0 = ynew;
                    n++;
                } while (x0 < x && n < short.MaxValue);
                result = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x0"></param>
        /// <param name="x"></param>
        /// <param name="h"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static double rungeKutta4(IDoubleMesh f, double x0, double x, double h, double y0)
        {
            int n = 0;
            double xnew, ynew, k1, k2, k3, k4, result = double.NaN;
            if (x == x0)
                result = y0;
            else if (x > x0)
            {
                do
                {
                    if (h > x - x0) h = x - x0;
                    k1 = h * f(x0, y0);
                    k2 = h * f(x0 + 0.5 * h, y0 + 0.5 * k1);
                    k3 = h * f(x0 + 0.5 * h, y0 + 0.5 * k2);
                    k4 = h * f(x0 + h, y0 + k3);
                    ynew = y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                    xnew = x0 + h;
                    x0 = xnew;
                    y0 = ynew;
                    n++;
                } while (x0 < x && n <  short.MaxValue);
                result = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x0"></param>
        /// <param name="x"></param>
        /// <param name="h"></param>
        /// <param name="y0"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        private static double rungeKuttaFehlberg(IDoubleMesh f, double x0, double x, double h, double y0, double tolerance)
        {
            int n = 0;
            double xnew, ynew, hnew, k1, k2, k3, k4, k5, k6;
            double hmin = 0.0001;
            double hmax = 0.5;
            if (h > hmax) h = hmax;
            if (h < hmin) h = hmin;

            while (x0 < x && n < short.MaxValue)
            {
                k1 = h * f(x0, y0);
                k2 = h * f(x0 + 0.25 * h, y0 + 0.25 * k1);
                k3 = h * f(x0 + 3 * h / 8, y0 + 3 * k1 / 32 + 9 * k2 / 32);
                k4 = h * f(x0 + 12 * h / 13, y0 + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197);
                k5 = h * f(x0 + h, y0 + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104);
                k6 = h * f(x0 + 0.5 * h, y0 - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40);
                double error = Math.Abs(k1 / 360 - 128 * k3 / 4275 - 2197 * k4 / 75240 + k5 / 50 + 2 * k6 / 55) / h;
                double s = Math.Pow(0.5 * tolerance / error, 0.25);
                if (error < tolerance)
                {
                    ynew = y0 + 25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - 0.2 * k5;
                    xnew = x0 + h;
                    x0 = xnew;
                    y0 = ynew;
                }
                if (s < 0.1) s = 0.1;
                if (s > 4) s = 4;
                hnew = h * s;
                h = hnew;
                if (h > hmax) h = hmax;
                if (h < hmin) h = hmin;
                if (h > x - x0) h = x - x0;
                n++;
            } return y0;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Метод дифферецирования.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Метод Эйлера.
            /// </summary>
            Euler,
            /// <summary>
            /// Метод Рунге-Кутты второго порядка.
            /// </summary>
            RungeKutta2,
            /// <summary>
            /// Метод Рунге-Кутты четвертого порядка.
            /// </summary>
            RungeKutta4,
            /// <summary>
            /// Метод Фелберга.
            /// </summary>
            Fehlberg,
            #endregion
        }
        #endregion
    }
    #endregion

    #region Approximation methods
    /// <summary>
    /// Определяет класс аппроксимации методом наименьших квадратов.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
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
        /// Инициализирует класс аппроксимации методом наименьших квадратов.
        /// </summary>
        /// <param name="power">Степень полинома</param>
        /// <param name="method">Метод аппроксимации</param>
        public Approximation(int power = 1, Approximation.Method method = Approximation.Method.Polynomial)
        {
            this.Power = power;
            this.method = method;
        }
        /// <summary>
        /// Получает или задает степень полинома.
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
                    throw new Exception("Неверное значение аргмуента");

                this.power = value;
            }
        }
        /// <summary>
        /// Получает или задает метод аппроксимации.
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
        /// Возвращает значение аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
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
        /// Возвращает значение аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="cf">Коэффициенты аппроксимации</param>
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
        /// Возвращает значение аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="cf">Коэффициенты аппроксимации</param>
        /// <param name="error">Погрешность аппроксимации</param>
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
        /// Возвращает значение аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="cf">Коэффициенты аппроксимации</param>
        /// <param name="error">Погрешность аппроксимации</param>
        /// <param name="equation">Уравнение аппроксимации</param>
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
        /// Возвращает значение аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
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
        /// Возвращает значение аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="cf">Коэффициенты аппроксимации</param>
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
        /// Возвращает значение аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="cf">Коэффициенты аппроксимации</param>
        /// <param name="error">Погрешность аппроксимации</param>
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
        /// Возвращает значение аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="cf">Коэффициенты аппроксимации</param>
        /// <param name="error">Погрешность аппроксимации</param>
        /// <param name="equation">Уравнение аппроксимации</param>
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
            cf = LeastSquares.Coefficients(x, y, m);
            double[] ya = LeastSquares.Polynomial(x, cf);
            error = LeastSquares.Error(ya, y);
            equation = LeastSquares.Equation(cf);
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
            cf = LeastSquares.Coefficients(x, y, m);
            Complex[] ya = LeastSquares.Polynomial(x, cf);
            error = LeastSquares.Error(ya, y);
            equation = LeastSquares.Equation(cf);
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
            cf = LeastSquares.Coefficients(xa, y, m);
            ya = LeastSquares.Polynomial(xa, cf);
            error = LeastSquares.Error(ya, y);
            equation = LeastSquares.Equation(cf, " * LN(X)^");
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
            cf = LeastSquares.Coefficients(xa, y, m);
            ya = LeastSquares.Polynomial(xa, cf);
            error = LeastSquares.Error(ya, y);
            equation = LeastSquares.Equation(cf, " * LN(X)^");
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
            cf = LeastSquares.Coefficients(x, ya, m);
            double[] p = LeastSquares.Polynomial(x, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Pow(Math.E, p[i]);
            }

            error = LeastSquares.Error(ya, y);
            equation = "EXP" + '(' + LeastSquares.Equation(cf) + ')';
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
            cf = LeastSquares.Coefficients(x, ya, m);
            Complex[] p = LeastSquares.Polynomial(x, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Pow(Math.E, p[i]);
            }

            error = LeastSquares.Error(ya, y);
            equation = "EXP" + '(' + LeastSquares.Equation(cf) + ')';
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
            cf = LeastSquares.Coefficients(xa, ya, m);
            double[] p = LeastSquares.Polynomial(xa, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Exp(p[i]);
            }

            error = LeastSquares.Error(ya, y);
            equation = "EXP" + '(' + LeastSquares.Equation(cf, " * LN(X)^") + ')';
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
            cf = LeastSquares.Coefficients(xa, ya, m);
            Complex[] p = LeastSquares.Polynomial(xa, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Maths.Exp(p[i]);
            }

            error = LeastSquares.Error(ya, y);
            equation = "EXP" + '(' + LeastSquares.Equation(cf, " * LN(X)^") + ')';
            return ya;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Метод аппроксимации.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Полиномиальная аппроксимация.
            /// </summary>
            Polynomial,
            /// <summary>
            /// Логарифимическая аппроксимация.
            /// </summary>
            Logarithmic,
            /// <summary>
            /// Экспоненциальная аппроксимация.
            /// </summary>
            Exponential,
            /// <summary>
            /// Степенная аппроксимация.
            /// </summary>
            Power,
            #endregion
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс, реализующий метод наименьших квадратов.
    /// </summary>
    internal static class LeastSquares
    {
        #region double components
        /// <summary>
        /// Возвращает значение полиномиала.
        /// </summary>
        /// <param name="x">Аргумент</param>
        /// <param name="c">Коэффициенты аппроксимации</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает массив значений полиномиала.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="c">Коэффициенты аппроксимации</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Polynomial(double[] x, double[] c)
        {
            int n = x.Length, i;
            double[] y = new double[n];

            for (i = 0; i < n; i++)
            {
                y[i] = LeastSquares.Polynomial(x[i], c);
            }
            return y;
        }
        /// <summary>
        /// Возвращает коэффициенты аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргмента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="iterations">Количество итераций</param>
        public static double[] Coefficients(double[] x, double[] y, int iterations)
        {
            // Построение матрицы преобразования:
            int i, j;
            int n = x.Length;
            int m = iterations < 1 ? 1 : iterations;
            double[,] matrix = new double[m, m + 1];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    matrix[i, j] = LeastSquares.SummaryPow(x, j + i);
                }
                matrix[i, m] = LeastSquares.SummaryPow(y, x, 1, i);
            }

            // Решение системы линейных уравнений:
            return Matrice.Solve(matrix);
        }
        /// <summary>
        /// Возвращает значение выражения: s += v(i) ^ pow.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="pow">Степень</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает значение выражения: s += {x(i) ^ powx} * {y(i) ^ powy}.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="y">Одномерный массив</param>
        /// <param name="powx">Степень</param>
        /// <param name="powy">Степень</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает погрешность аппроксимации функции.
        /// </summary>
        /// <param name="a">Аппроксимация</param>
        /// <param name="b">Функция</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает уравнение полинома, представленного в виде строки.
        /// </summary>
        /// <param name="p">Коэффициенты полинома</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
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
        /// Возвращает уравнение полинома, представленного в виде строки.
        /// </summary>
        /// <param name="p">Коэффициенты полинома</param>
        /// <param name="function">Функция</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
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
        /// Возвращает значение полиномиала.
        /// </summary>
        /// <param name="x">Аргумент</param>
        /// <param name="c">Коэффициенты аппроксимации</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает массив значений полиномиала.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="c">Коэффициенты аппроксимации</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Polynomial(Complex[] x, Complex[] c)
        {
            int n = x.Length, i;
            Complex[] y = new Complex[n];

            for (i = 0; i < n; i++)
            {
                y[i] = LeastSquares.Polynomial(x[i], c);
            }
            return y;
        }
        /// <summary>
        /// Возвращает коэффициенты аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргмента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="iterations">Количество итераций</param>
        public static Complex[] Coefficients(Complex[] x, Complex[] y, int iterations)
        {
            // Построение матрицы преобразования:
            int i, j;
            int n = x.Length;
            int m = iterations < 1 ? 1 : iterations;
            Complex[,] matrix = new Complex[m, m + 1];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    matrix[i, j] = LeastSquares.SummaryPow(x, j + i);
                }
                matrix[i, m] = LeastSquares.SummaryPow(y, x, 1, i);
            }

            // Решение системы линейных уравнений:
            return Matrice.Solve(matrix);
        }
        /// <summary>
        /// Возвращает значение выражения: s += v(i) ^ pow.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="pow">Степень</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает значение выражения: s += {x(i) ^ powx} * {y(i) ^ powy}.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="y">Одномерный массив</param>
        /// <param name="powx">Степень</param>
        /// <param name="powy">Степень</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает погрешность аппроксимации функции.
        /// </summary>
        /// <param name="a">Аппроксимация</param>
        /// <param name="b">Функция</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает уравнение полинома, представленного в виде строки.
        /// </summary>
        /// <param name="p">Коэффициенты полинома</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
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
        /// Возвращает уравнение полинома, представленного в виде строки.
        /// </summary>
        /// <param name="p">Коэффициенты полинома</param>
        /// <param name="function">Функция</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
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
    /// Определяет класс решения уравнений с использованием спектрального разложения матрицы.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://www.mathworks.com/help/matlab/ref/roots.html
    /// </remarks>
    /// </summary>
    public class Roots
    {
        #region Private data
        /// <summary>
        /// Спектральное разложение матрицы.
        /// </summary>
        private EVD eig;
        /// <summary>
        /// Погрешность вычислений.
        /// </summary>
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс решения уравнений с использованием спектрального разложения матрицы.
        /// </summary>
        /// <param name="eps">Погрешность [0, 1]</param>
        public Roots(double eps = 1e-16)
        {
            this.Eps = eps;
        }
        /// <summary>
        /// Получает или задает погрешность [0, 1].
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
        /// Возвращает вектор-столбец, соотвествующий численному решению полинома: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="polynomial">Полином</param>
        /// <returns>Вектор-столбец</returns>
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
        /// Возвращает вектор-столбец коэффициентов полинома: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="roots">Корни полинома</param>
        /// <returns>Вектор-столбец</returns>
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
