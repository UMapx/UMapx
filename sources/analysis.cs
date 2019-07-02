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
    //                               UMAPX.ANALYSIS
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Integral solution
    /// <summary>
    /// Определяет класс, реализующий численное интегрирование методом трапеций.
    /// </summary>
    public class TrapezoidIntegralSolution : INumericSolution
    {
        #region Private data
        /// <summary>
        /// Делегат непрерывной функции.
        /// </summary>
        private IDouble f;
        /// <summary>
        /// Погрешность вычислений.
        /// </summary>
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий численное интегрирование методом трапеций.
        /// </summary>
        public TrapezoidIntegralSolution() { }
        /// <summary>
        /// Инициализирует класс, реализующий численное интегрирование методом трапеций.
        /// </summary>
        /// <param name="function">Делегат непрерывной фукнции</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public TrapezoidIntegralSolution(IDouble function, double eps = 0.0001)
        {
            this.f = function;
            this.Eps = eps;
        }
        /// <summary>
        /// Получает или задает делегат непрерывной функции.
        /// </summary>
        public IDouble Function
        {
            get
            {
                return this.f;
            }
            set
            {
                this.f = value;
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
        /// Возвращает значение интеграла функции.
        /// </summary>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(double a, double b)
        {
            // Trapezoid Integral Solution:
            // by Asiryan Valeriy

            // compute first points:
            int n = 1, i;
            double h = (b - a);
            double s = (f(a) + f(b)) * h / 2.0;
            double s1 = 0;

            // process:
            while (Math.Abs(s - s1) > eps)
            {
                n = n * 2;
                h = (b - a) / n;
                s1 = s;
                s = f(a) + f(b);

                for (i = 0; i < n; i++)
                {
                    s += 2 * f(a + i * h);
                }
                s *= (h / 2.0);
            }

            return s;
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс, реализующий численное интегрирование методом Симпсона.
    /// </summary>
    public class SimpsonIntegralSolution : INumericSolution
    {
        #region Private data
        /// <summary>
        /// Делегат непрерывной функции.
        /// </summary>
        private IDouble f;
        /// <summary>
        /// Погрешность вычислений.
        /// </summary>
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий численное интегрирование методом Симпсона.
        /// </summary>
        public SimpsonIntegralSolution() { }
        /// <summary>
        /// Инициализирует класс, реализующий численное интегрирование методом Симпсона.
        /// </summary>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public SimpsonIntegralSolution(IDouble function, double eps = 0.0001) 
        { 
            this.f = function;
            this.Eps = eps;
        }
        /// <summary>
        /// Получает или задает делегат непрерывной функции.
        /// </summary>
        public IDouble Function
        {
            get
            {
                return this.f;
            }
            set
            {
                this.f = value;
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
        /// Возвращает значение интеграла функции.
        /// </summary>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(double a, double b)
        {
            // Simpson Integral Solution:
            // by Angelina Drabkova
            int c, n, k;
            double s1 = 0, x, s, h;

            // compute first points:
            n = 2; h = (b - a) / n;
            s = s = (f(a) + 4 * f((a + b) / 2) + f(b)) * h / 3;

            // process:
            while (Maths.Abs(s - s1) > eps)
            {
                h = h / 2.0;
                n = 2 * n;
                s1 = s;
                c = 4;
                x = a;
                s = f(a) + f(b);

                for (k = 0; k < n - 1; k++)
                {
                    x = x + h;
                    s = s + c * f(x);
                    c = 6 - c;
                }
                s = s * h / 3.0;
            }
            return s;
        }
        #endregion
    }
    #endregion

    #region Differential equation solution
    /// <summary>
    /// Определяет класс, реализующий решение дифференциального уравнения методом Рунге-Кутты.
    /// </summary>
    public class RungeKuttaSolution : INumericSolution
    {
        #region Private data
        /// <summary>
        /// Делегат непрерывной функции.
        /// </summary>
        private IDoubleMesh f;
        /// <summary>
        /// Погрешность вычислений.
        /// </summary>
        private double eps;
        /// <summary>
        /// Значение функции.
        /// </summary>
        private double y0;
        #endregion

        #region Runge-Kutta components
        /// <summary>
        /// Инициализирует класс, реализующий решение дифференциального уравнения методом Рунге-Кутты.
        /// </summary>
        /// <param name="function">Делегат непрерывной функции, зависящей от двух переменных</param>
        /// <param name="y0">Значение функции</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public RungeKuttaSolution(IDoubleMesh function, double y0, double eps = 0.0001)
        {
            Function = function;
            Y0 = y0;
            Eps = eps;
        }
        /// <summary>
        /// Получает или задает делегат непрерывной функции, зависящей от двух переменных.
        /// </summary>
        public IDoubleMesh Function
        {
            get
            {
                return this.f;
            }
            set
            {
                this.f = value;
            }
        }
        /// <summary>
        /// Получает или задает значение функции.
        /// </summary>
        public double Y0
        {
            get
            {
                return this.y0;
            }
            set
            {
                this.y0 = value;
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
        /// Возвращает значение дифференциального уравнения.
        /// </summary>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(double a, double b)
        {
            // Meain Runge-Kutta algorithm
            // by T.I. Semenova:

            int m = 1;
            double h = b - a;
            double y = rk4(a, y0, h, m);
            double z = double.MaxValue;

            while (Math.Abs(y - z) > eps)
            {
                z = y;
                y = y0;
                h = h / 2;
                m = 2 * m;
                y = rk4(a, y0, h, m);
            }

            return y;
        }
        /// <summary>
        /// Реализует алгоритм Рунге-Кутты для решения дифференциального уравнения.
        /// </summary>
        /// <param name="x0">Аргумент</param>
        /// <param name="y0">Функция</param>
        /// <param name="h">Шаг</param>
        /// <param name="m">Порядок</param>
        /// <returns>Значение функции</returns>
        private double rk4(double x0, double y0, double h, int m)
        {
            // Runge-Kutta iterative algorithm
            // by T.I. Semenova:

            double x = x0, y = y0;
            double k1, k2, k3, k4;

            for (int j = 0; j < m; j++)
            {
                k1 = f(x, y);
                k2 = f(x + h / 2, y + h * k1 / 2);
                k3 = f(x + h / 2, y + h * k2 / 2);
                k4 = f(x + h, y + h * k3);
                y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                x = x + h;
            }

            return y;
        }
        #endregion
    }
    #endregion

    #region Nonlinear solution
    /// <summary>
    /// Определяет класс, реализующий решение нелинейного уравнения методом половинного деления.
    /// </summary>
    public class BisectionNonlinearSolution : INumericSolution
    {
        #region Private data
        /// <summary>
        /// Делегат непрерывной функции.
        /// </summary>
        private IDouble f;
        /// <summary>
        /// Погрешность вычислений.
        /// </summary>
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий решение нелинейного уравнения методом половинного деления.
        /// </summary>
        public BisectionNonlinearSolution() { }
        /// <summary>
        /// Инициализирует класс, реализующий решение нелинейного уравнения методом половинного деления.
        /// </summary>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public BisectionNonlinearSolution(IDouble function, double eps = 0.0001)
        {
            this.f = function;
            this.Eps = eps;
        }
        /// <summary>
        /// Получает или задает делегат непрерывной функции.
        /// </summary>
        public IDouble Function
        {
            get
            {
                return this.f;
            }
            set
            {
                this.f = value;
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
        /// <param name="a">Начало отрезка</param>
        /// <param name="b">Конец отрезка</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(double a, double b)
        {
            double c = (b + a) / 2.0;
            int n = 0;

            while ((b - a) > eps)
            {
                n++;

                if ((f(c) * f(b)) <= 0)
                {
                    a = c;
                }
                else
                {
                    b = c;
                }
                c = (a + b) / 2;
            }
            return c;
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс, реализующий решение нелинейного уравнения методом хорд.
    /// </summary>
    public class ChordNonlinearSolution : INumericSolution
    {
        #region Private data
        /// <summary>
        /// Делегат непрерывной функции.
        /// </summary>
        private IDouble f;
        /// <summary>
        /// Погрешность вычислений.
        /// </summary>
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий решение нелинейного уравнения методом хорд.
        /// </summary>
        public ChordNonlinearSolution() { }
        /// <summary>
        /// Инициализирует класс, реализующий решение нелинейного уравнения методом хорд.
        /// </summary>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public ChordNonlinearSolution(IDouble function, double eps = 0.0001)
        {
            this.f = function;
            this.Eps = eps;
        }
        /// <summary>
        /// Получает или задает делегат непрерывной функции.
        /// </summary>
        public IDouble Function
        {
            get
            {
                return this.f;
            }
            set
            {
                this.f = value;
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
        /// Возвращает значения корня нелинейного уравнения.
        /// </summary>
        /// <param name="a">Начало отрезка</param>
        /// <param name="b">Конец отрезка</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(double a, double b)
        {
            double x0 = (b - a) / 2.0;
            double x;

            while (Maths.Abs(f(x0) / b) > eps)
            {
                x = x0;
                x0 = x - (f(x) * (a - x)) / (f(a) - f(x));
            }
            return x0;
        }
        #endregion
    }
    #endregion

    #region Interpolation methods
    /// <summary>
    /// Определяет класс для интерполяции методом Лагранжа.
    /// </summary>
    public class LagrangeInterpolation : IInterpolation
    {
        #region Private data
        /// <summary>
        /// Массив значений аргумента.
        /// </summary>
        private double[] x;
        /// <summary>
        /// Массив значений функции.
        /// </summary>
        private double[] y;
        /// <summary>
        /// Количество итераций.
        /// </summary>
        private int iterations;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс для интерполяции методом Лагранжа.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="iterations">Количество итераций</param>
        public LagrangeInterpolation(double[] x, double[] y, int iterations)
        {
            X = x; Y = y; Iterations = iterations;
        }
        /// <summary>
        /// Получает или задает массив значений аргумента.
        /// </summary>
        public double[] X
        {
            get
            {
                return this.x;
            }
            set
            {
                this.x = value;
            }
        }
        /// <summary>
        /// Получает или задает массив значений функции.
        /// </summary>
        public double[] Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = value;
            }
        }
        /// <summary>
        /// Получает или задает количество итераций.
        /// </summary>
        public int Iterations
        {
            get
            {
                return this.iterations;
            }
            set
            {
                this.iterations = value;
            }
        }
        /// <summary>
        /// Возвращает значение функции в точке.
        /// </summary>
        /// <param name="xl">Значение аргумента для вычисления</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(double xl)
        {
            // Lagrange interpolation algorithm
            // by Asiryan Valeriy.

            int k = Maths.Range(iterations, 1, y.Length - 1);
            double l = 0, l1;
            int i, j;

            for (i = 0; i <= k; i++)
            {
                l1 = 1;

                for (j = 0; j <= k; j++)
                {
                    if (j != i)
                    {
                        l1 = (xl - x[j]) / (x[i] - x[j]) * l1;
                    }
                }
                l = l + l1 * y[i];
            }

            return l;
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс для интерполяции по первой формуле Ньютона.
    /// </summary>
    public class NewtonInterpolation : IInterpolation
    {
        #region Private data
        /// <summary>
        /// Массив значений аргумента.
        /// </summary>
        private double[] x;
        /// <summary>
        /// Массив значений функции.
        /// </summary>
        private double[] y;
        /// <summary>
        /// Количество итераций.
        /// </summary>
        private int iterations;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс для интерполяции по первой формуле Ньютона.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="iterations">Количество итераций</param>
        public NewtonInterpolation(double[] x, double[] y, int iterations = 10)
        {
            X = x; Y = y; Iterations = iterations;
        }
        /// <summary>
        /// Получает или задает массив значений аргумента.
        /// </summary>
        public double[] X
        {
            get
            {
                return this.x;
            }
            set
            {
                this.x = value;
            }
        }
        /// <summary>
        /// Получает или задает массив значений функции.
        /// </summary>
        public double[] Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = value;
            }
        }
        /// <summary>
        /// Получает или задает количество итераций.
        /// </summary>
        public int Iterations
        {
            get
            {
                return this.iterations;
            }
            set
            {
                this.iterations = value;
            }
        }
        /// <summary>
        /// Возвращает значение функции в точке.
        /// </summary>
        /// <param name="xl">Значение аргумента для вычисления</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(double xl)
        {
            // Newton interpolation algorithm
            // by Asiryan Valeriy.

            int n = Maths.Range(iterations, 1, x.Length);
            double res = y[0], F, den;
            int i, j, k;

            for (i = 1; i < n; i++)
            {
                F = 0;
                for (j = 0; j <= i; j++)
                {
                    den = 1;
                    for (k = 0; k <= i; k++)
                    {
                        if (k != j) den *= (x[j] - x[k]);
                    }
                    F += y[j] / den;
                }
                for (k = 0; k < i; k++) F *= (xl - x[k]);
                res += F;
            }
            return res;
        }
        #endregion
    }
    #endregion

    #region Optimization methods
    /// <summary>
    /// Определяет класс, реализующий поиск минимума методом золотого сечения.
    /// </summary>
    public class GoldenOptimization : INumericSolution
    {
        #region Private data
        /// <summary>
        /// Делегат непрерывной функции.
        /// </summary>
        private IDouble f;
        /// <summary>
        /// Погрешность вычислений.
        /// </summary>
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий поиск минимума методом золотого сечения.
        /// </summary>
        public GoldenOptimization() { }
        /// <summary>
        /// Инициализирует класс, реализующий поиск минимума методом золотого сечения.
        /// </summary>
        /// <param name="function">Делегат непрерывной фукнции</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public GoldenOptimization(IDouble function, double eps = 0.0001)
        {
            this.f = function;
            this.Eps = eps;
        }
        /// <summary>
        /// Получает или задает делегат непрерывной функции.
        /// </summary>
        public IDouble Function
        {
            get
            {
                return this.f;
            }
            set
            {
                this.f = value;
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
        /// Возвращает соответствующий минимум функции на отреке.
        /// </summary>
        /// <param name="a">Начало отрезка</param>
        /// <param name="b">Конец отрезка</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Compute(double a, double b)
        {
            double lbound = a, rbound = b;
            double l = rbound - (rbound - lbound) / Maths.Phi;
            double r = lbound + (rbound - lbound) / Maths.Phi;

            while (Maths.Abs(rbound - lbound) >= eps)
            {
                if (f(l) >= f(r))
                {
                    lbound = l;
                    l = r;
                    r = lbound + (rbound - lbound) / Maths.Phi;
                }
                else
                {
                    rbound = r;
                    r = l;
                    l = rbound - (rbound - lbound) / Maths.Phi;
                }
            }

            return (lbound + rbound) / 2.0;
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс, реализующий поиск минимума методом дихотомии.
    /// </summary>
    public class DichotomyOptimization : INumericSolution
    {
        #region Private data
        /// <summary>
        /// Делегат непрерывной функции.
        /// </summary>
        private IDouble f;
        /// <summary>
        /// Погрешность вычислений.
        /// </summary>
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий поиск минимума методом дихотомии.
        /// </summary>
        public DichotomyOptimization() { }
        /// <summary>
        /// Инициализирует класс, реализующий поиск минимума методом дихотомии.
        /// </summary>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public DichotomyOptimization(IDouble function, double eps = 0.0001)
        {
            this.f = function;
            this.Eps = eps;
        }
        /// <summary>
        /// Получает или задает делегат непрерывной функции.
        /// </summary>
        public IDouble Function
        {
            get
            {
                return this.f;
            }
            set
            {
                this.f = value;
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
        /// Возвращает соответствующий минимум функции на отреке.
        /// </summary>
        /// <param name="c">Нижний предел</param>
        /// <param name="d">Верхний предел</param>
        /// <returns>Пара дробных чисел X, Y</returns>
        public double Compute(double c, double d)
        {
            double a = c;
            double b = d;
            double del, x1, x2;
            del = eps / 10.0;

            while (Math.Abs(b - a) >= eps)
            {
                x1 = (a + b - del) / 2.0;
                x2 = (a + b + del) / 2.0;

                if (f(x1) > f(x2))
                {
                    a = x1;
                }
                else
                {
                    b = x2;
                }
            }

            return (a + b) / 2.0;
        }
        #endregion
    }
    #endregion

    #region SLAU methods
    /// <summary>
    /// Определяет класс, реализующий решение системы линейных алгебраических уравнений методом Якоби.
    /// <remarks>
    /// В числовой линейной алгебре метод Якоби (или итерационный метод Якоби) является алгоритмом определения решений диагонально доминирующей системы
    /// линейных уравнений. Процесс вычислений повторяется до тех пор, пока он предел не сходится. Этот алгоритм является упрощенной версией метода 
    /// преобразования Якоби матричной диагонализации.
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Jacobi_method
    /// </remarks>
    /// </summary>
    public class JacobiSolution : IMatrixSolution
    {
        #region Private data
        private double[] b;
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий решение системы линейных алгебраических уравнений методом Якоби.
        /// </summary>
        /// <param name="b">Вектор-строка</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public JacobiSolution(double[] b, double eps = 0.0001)
        {
            B = b; Eps = eps;
        }
        /// <summary>
        /// Получает или задает вектор-строку.
        /// </summary>
        public double[] B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
            }
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
        /// Возвращает вектор-строку, соотвествующий решению системы линейных алгебраических уравнений: Ax = b.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        /// <returns>Вектор-строка</returns>
        public double[] Compute(double[,] A)
        {
            int N = b.Length;
            double[] p = new double[N];
            double[] x = new double[N];
            double norm;
            int i, j, k;

            do
            {
                for (i = 0; i < N; i++)
                {
                    p[i] = b[i];
                    for (j = 0; j < N; j++)
                    {
                        if (i != j)
                        {
                            p[i] -= A[i, j] * x[j];
                        }
                    }
                    p[i] /= A[i, i];
                }
                norm = Math.Abs(x[0] - p[0]);
                for (k = 0; k < N; k++)
                {
                    if (Math.Abs(x[k] - p[k]) > norm)
                    {
                        norm = Math.Abs(x[k] - p[k]);
                    }
                    x[k] = p[k];
                }
            } while (norm > eps);

            return x;
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс, реализующий решение системы линейных алгебраических уравнений методом Гаусса-Зейделя.
    /// <remarks>
    /// В числовой линейной алгебре метод Гаусса-Зейделя, также известный как метод Либмана или метод последовательного перемещения, 
    /// является итерационным методом, используемым для решения линейной системы уравнений. Он назван в честь немецких математиков Карла Фридриха Гаусса
    /// и Филиппа Людвига фон Зейделя и похож на метод Якоби. Хотя его можно применять к любой матрице с ненулевыми элементами на диагоналях, сходимость 
    /// гарантируется только в том случае, если матрица либо диагонально доминантна, либо симметрична и положительно определена.
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
    /// </remarks>
    /// </summary>
    public class GaussSeidelSolution : IMatrixSolution
    {
        #region Private data
        private double[] b;
        private double eps;
        #endregion
        
        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий решение системы линейных алгебраических уравнений методом Гаусса-Зейделя.
        /// </summary>
        /// <param name="b">Вектор-строка</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public GaussSeidelSolution(double[] b, double eps = 0.0001)
        {
            B = b; Eps = eps;
        }
        /// <summary>
        /// Получает или задает вектор-строку.
        /// </summary>
        public double[] B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Получает или задает погрешность [0, 1]
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
        /// Возвращает вектор-строку, соотвествующий решению системы линейных алгебраических уравнений: Ax = b.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        /// <returns>Вектор-строка</returns>
        public double[] Compute(double[,] A)
        {
            int n = b.Length;
            double[] x = new double[n];
            double[] temp;
            int i, j;
            double var;

            do
            {
                temp = (double[])x.Clone();

                for (i = 0; i < n; i++)
                {
                    var = 0;
                    for (j = 0; j < i; j++)
                    {
                        var += (A[i, j] * x[j]);
                    }
                    for (j = i + 1; j < n; j++)
                    {
                        var += (A[i, j] * temp[j]);
                    }
                    x[i] = (b[i] - var) / A[i, i];
                }
            } while (!Converge(x, temp, eps, n));

            return x;
        }
        /// <summary>
        /// Возвращает логическое значение при условии сходимости.
        /// </summary>
        /// <param name="xk">Вектор-строка</param>
        /// <param name="xkp">Вектор-строка</param>
        /// <param name="eps">Погрешность вычислений</param>
        /// <param name="n">Длина вектора</param>
        /// <returns>Логическое значение</returns>
        private bool Converge(double[] xk, double[] xkp, double eps, int n)
        {
            double norm = 0;
            for (int i = 0; i < n; i++)
            {
                norm += (xk[i] - xkp[i]) * (xk[i] - xkp[i]);
            }
            if (Math.Sqrt(norm) >= eps)
            {
                return false;
            }
            return true;
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс, реализующий решение квадратных систем линейных алгебраических уравнений методом Гаусса-Йордана (метод полного исключения неизвестных).
    /// <remarks>
    /// Метод является модификацией метода Гаусса. Назван в честь Карла Фридриха Гаусса и Вильгельма Йордана. 
    /// Этот метод также может быть использован для нахождения ранга матрицы, вычисления определителя матрицы и вычисления обратного к обратимой квадратной матрице. 
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Gauss%E2%80%93Jordan_elimination
    /// </remarks>
    /// </summary>
    public class GaussJordanElimination : IMatrixSolution
    {
        #region Class components
        /// <summary>
        /// Инициализирует класс, реализующий решение квадратных систем линейных алгебраических уравнений методом Гаусса-Йордана.
        /// </summary>
        public GaussJordanElimination() { }
        /// <summary>
        /// Возвращает вектор-строку, соотвествующий решению системы линейных алгебраических уравнений: Ax = b.
        /// <remarks>
        /// По завершению процедуры исходная расширенная матрица A будет представлять собой верхнеугольную матрицу.
        /// </remarks>
        /// </summary>
        /// <param name="A">Расширенная матрица</param>
        /// <returns>Вектор-строка</returns>
        public double[] Compute(double[,] A)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            int i, j, k, wf, hf, l = width - 1;
            double[] x = new double[height];
            double temp;

            // Этап 1:
            for (i = 0; i < height; i++)
            {
                // Делаем главную диагональ единицами:
                temp = A[i, i];
                for (j = 0; j < width; j++)
                {
                    A[i, j] /= temp;
                }
                // Обнуляем числа под единицами главной диогoнали:
                for (k = i + 1; k < height; k++)
                {
                    temp = A[k, i];
                    for (j = i; j < width; j++)
                    {
                        A[k, j] = A[k, j] - A[i, j] * temp;
                    }
                }
            }

            // Этап 2:
            for (i = 0; i < height; i++)
            {
                wf = (width - 1) - i - 1;
                hf = (height - 1) - i;

                for (k = 0; k < (height - 1) - i; k++)
                {
                    temp = A[k, wf];

                    for (j = wf; j < width; j++)
                    {
                        A[k, j] = A[k, j] - A[hf, j] * temp;
                    }
                }
            }

            for (k = 0; k < height; k++)
            {
                x[k] = A[k, l];
            }
            return x;
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
    public class EigenRootsSolution : IRootsSolution
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
        public EigenRootsSolution(double eps = 1e-16)
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
        /// Возвращает вектор-строку, соотвествующий численному решению полинома: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="polynomial">Полином</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор-строку коэффициентов полинома: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="roots">Корни полинома</param>
        /// <returns>Одномерный массив</returns>
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

    #region Approximation methods
    /// <summary>
    /// Определяет класс полиномиальной аппроксимации.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// http://simenergy.ru/math-analysis/digital-processing/85-ordinary_least_squares
    /// </remarks>
    /// </summary>
    public class PolynomialApproximation : IApproximation
    {
        #region Private data
        double[] x, y, ya, cf;
        #endregion

        #region Approximation components
        /// <summary>
        /// Инициализирует класс полиномиальной аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="power">Степень полинома</param>
        public PolynomialApproximation(double[] x, double[] y, int power = 1)
        {
            // Options:
            this.x = x; this.y = y;
            int m = (power < 1) ? 2 : power + 1;
            this.cf = LeastSquares.Coefficients(x, y, m);
            this.ya = LeastSquares.Polynomial(x, cf);
            return;
        }
        /// <summary>
        /// Возвращает значение погрешности аппроксимации.
        /// </summary>
        public double Error
        {
            get
            {
                return LeastSquares.Error(ya, y);
            }
        }
        /// <summary>
        /// Возвращает коэффициенты полинома.
        /// </summary>
        public double[] Coefficients
        {
            get
            {
                return this.cf;
            }
        }
        /// <summary>
        /// Возвращает аппроксимацию функции.
        /// </summary>
        public double[] Approximation
        {
            get
            {
                return this.ya;
            }
        }
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает уравнение аппроксимации.
        /// </summary>
        /// <returns>Текст как последовательность знаков юникода</returns>
        public override string ToString()
        {
            return LeastSquares.Equation(this.cf);
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс логарифмической аппроксимации.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// http://simenergy.ru/math-analysis/digital-processing/85-ordinary_least_squares
    /// </remarks>
    /// </summary>
    public class LogarithmicApproximation : IApproximation
    {
        #region Private data
        double[] x, y, xa, ya, cf;
        #endregion

        #region Approximation components
        /// <summary>
        /// Инициализирует класс логарифмической аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="power">Степень полинома</param>
        public LogarithmicApproximation(double[] x, double[] y, int power = 1)
        {
            // Options:
            this.x = x; this.y = y;
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            this.xa = new double[n];
            this.ya = new double[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Math.Log(x[i]);
            }

            // approximation:
            this.cf = LeastSquares.Coefficients(xa, y, m);
            this.ya = LeastSquares.Polynomial(xa, cf);
            return;
        }
        /// <summary>
        /// Возвращает значение погрешности аппроксимации.
        /// </summary>
        public double Error
        {
            get
            {
                return LeastSquares.Error(ya, y);
            }
        }
        /// <summary>
        /// Возвращает коэффициенты полинома.
        /// </summary>
        public double[] Coefficients
        {
            get
            {
                return this.cf;
            }
        }
        /// <summary>
        /// Возвращает аппроксимацию функции.
        /// </summary>
        public double[] Approximation
        {
            get
            {
                return this.ya;
            }
        }
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает уравнение аппроксимации.
        /// </summary>
        /// <returns>Текст как последовательность знаков юникода</returns>
        public override string ToString()
        {
            return LeastSquares.Equation(this.cf, " * LN(X)^");
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс экспоненциальной аппроксимации.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// http://simenergy.ru/math-analysis/digital-processing/85-ordinary_least_squares
    /// </remarks>
    /// </summary>
    public class ExponentialApproximation : IApproximation
    {
        #region Private data
        double[] x, y, ya, cf;
        #endregion

        #region Approximation components
        /// <summary>
        /// Инициализирует класс экспоненциальной аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="power">Степень полинома</param>
        public ExponentialApproximation(double[] x, double[] y, int power = 1)
        {
            // Options:
            this.x = x; this.y = y;
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            this.ya = new double[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Math.Log(y[i], Math.E);
            }

            // approximation:
            this.cf = LeastSquares.Coefficients(x, ya, m);
            double[] p = LeastSquares.Polynomial(x, this.cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Math.Pow(Math.E, p[i]);
            }

            return;
        }
        /// <summary>
        /// Возвращает значение погрешности аппроксимации.
        /// </summary>
        public double Error
        {
            get
            {
                return LeastSquares.Error(ya, y);
            }
        }
        /// <summary>
        /// Возвращает коэффициенты полинома.
        /// </summary>
        public double[] Coefficients
        {
            get
            {
                return this.cf;
            }
        }
        /// <summary>
        /// Возвращает аппроксимацию функции.
        /// </summary>
        public double[] Approximation
        {
            get
            {
                return this.ya;
            }
        }
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает уравнение аппроксимации.
        /// </summary>
        /// <returns>Текст как последовательность знаков юникода</returns>
        public override string ToString()
        {
            return "EXP" + '(' + LeastSquares.Equation(this.cf) + ')';
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс степенной аппроксимации.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// http://simenergy.ru/math-analysis/digital-processing/85-ordinary_least_squares
    /// </remarks>
    /// </summary>
    public class PowerApproximation : IApproximation
    {
        #region Private data
        double[] x, y, xa, ya, cf;
        #endregion

        #region Approximation components
        /// <summary>
        /// Инициализирует класс степенной аппроксимации.
        /// </summary>
        /// <param name="x">Массив значений аргумента</param>
        /// <param name="y">Массив значений функции</param>
        /// <param name="power">Степень полинома</param>
        public PowerApproximation(double[] x, double[] y, int power = 1)
        {
            // Options:
            this.x = x; this.y = y;
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            this.xa = new double[n];
            this.ya = new double[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = Math.Log(x[i]);
                ya[i] = Math.Log(y[i]);
            }

            // approximation:
            this.cf = LeastSquares.Coefficients(xa, ya, m);
            double[] p = LeastSquares.Polynomial(xa, this.cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = Math.Exp(p[i]);
            }
            return;
        }
        /// <summary>
        /// Возвращает значение погрешности аппроксимации.
        /// </summary>
        public double Error
        {
            get
            {
                return LeastSquares.Error(ya, y);
            }
        }
        /// <summary>
        /// Возвращает коэффициенты полинома.
        /// </summary>
        public double[] Coefficients
        {
            get
            {
                return this.cf;
            }
        }
        /// <summary>
        /// Возвращает аппроксимацию функции.
        /// </summary>
        public double[] Approximation
        {
            get
            {
                return this.ya;
            }
        }
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает уравнение аппроксимации.
        /// </summary>
        /// <returns>Текст как последовательность знаков юникода</returns>
        public override string ToString()
        {
            return "EXP" + '(' + LeastSquares.Equation(this.cf, " * LN(X)^") + ')';
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс, реализующий метод наименьших квадратов.
    /// </summary>
    internal static class LeastSquares
    {
        #region Private data
        /// <summary>
        /// Метод Гаусса-Йордана.
        /// </summary>
        private static GaussJordanElimination gje = new GaussJordanElimination();
        #endregion

        #region Static components
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
            return gje.Compute(matrix);
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
            if (x.Length != y.Length)
                throw new Exception("Размерности векторов должны быть одинаковыми");

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
    }
    #endregion

    #region Interfaces
    /// <summary>
    /// Определяет общий интерфейс численных методов.
    /// </summary>
    public interface INumericSolution
    {
        #region Components
        /// <summary>
        /// Возвращает значение функции.
        /// </summary>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        double Compute(double a, double b);
        #endregion
    }
    /// <summary>
    /// Определяет общий интерфейс решения уравнений.
    /// </summary>
    public interface IRootsSolution
    {
        #region Components
        /// <summary>
        /// Возвращает вектор-строку, соотвествующий численному полинома: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="polynomial">Полином</param>
        /// <returns>Одномерный массив</returns>
        Complex[] Compute(double[] polynomial);
        /// <summary>
        /// Возвращает вектор-строку коэффициентов полинома: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="roots">Корни полинома</param>
        /// <returns>Одномерный массив</returns>
        double[] Compute(Complex[] roots);
        #endregion
    }
    /// <summary>
    /// Определяет общий интерфейс решений системы линейных алгебраических уравнений.
    /// </summary>
    public interface IMatrixSolution
    {
        #region Components
        /// <summary>
        /// Возвращает вектор-строку, соотвествующий решению системы линейных алгебраических уравнений: Ax = b.
        /// </summary>
        /// <param name="A">Матрица</param>
        /// <returns>Вектор-строка</returns>
        double[] Compute(double[,] A);
        #endregion
    }
    /// <summary>
    /// Определяет общий интерфейс методов интерполяции.
    /// </summary>
    public interface IInterpolation
    {
        #region Components
        /// <summary>
        /// Получает или задает массив значений аргумента.
        /// </summary>
        double[] X { get; set; }
        /// <summary>
        /// Получает или задает массив значений функции.
        /// </summary>
        double[] Y { get; set; }
        /// <summary>
        /// Получает или задает количество итераций.
        /// </summary>
        int Iterations { get; set; }
        /// <summary>
        /// Возвращает значение функции в точке.
        /// </summary>
        /// <param name="xl">Значение аргумента для вычисления</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        double Compute(double xl);
        #endregion
    }
    /// <summary>
    /// Определяет общий интерфейс методов аппроксимации.
    /// </summary>
    public interface IApproximation
    {
        #region Components
        /// <summary>
        /// Возвращает значение погрешности аппроксимации.
        /// </summary>
        double Error { get; }
        /// <summary>
        /// Возвращает коэффициенты полинома.
        /// </summary>
        double[] Coefficients { get; }
        /// <summary>
        /// Возвращает аппроксимацию функции.
        /// </summary>
        double[] Approximation { get; }
        #endregion
    }
    #endregion
}
