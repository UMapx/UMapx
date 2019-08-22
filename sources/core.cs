// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Collections.Generic;
using System.Linq;

namespace UMapx.Core
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                                 UMAPX.CORE
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Mathematics modules
    /// <summary>
    /// Используется для реализации основных алгебраических, тригонометрических и гиперболических операций.
    /// </summary>
    public static class Maths
    {
        #region Constant
        /// <summary>
        /// Экспонента.
        /// </summary>
        public const double E = 2.7182818284590452353602874713527;
        /// <summary>
        /// Пи.
        /// </summary>
        public const double Pi = 3.141592653589793238462643383279;
        /// <summary>
        /// Фи (золотое число).
        /// </summary>
        public const double Phi = 1.6180339887498948482;
        /// <summary>
        /// Два пи.
        /// </summary>
        public const double Tau = 6.283185307179586476925286766558;
        /// <summary>
        /// Постоянная Эйлера-Маскерони.
        /// </summary>
        public const double Gamma = 0.577215664901532860606512090;
        /// <summary>
        /// Кадратный корень из числа 2.
        /// </summary>
        public const double Sqrt2 = 1.4142135623730950488016887242097;
        /// <summary>
        /// Постоянная Каталана.
        /// </summary>
        public const double G = 0.915965594177219015054603514932384110774;
        /// <summary>
        /// Постоянная Апери.
        /// </summary>
        public const double A = 1.202056903159594285399738161511449990764;
        /// <summary>   
        /// Мнимая единица.
        /// </summary>
        public static readonly Complex I = Complex.I;
        #endregion

        #region Types and ranges
        /// <summary>
        /// Преобразовывает случайную величину в тип Byte.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Целое число без знака</returns>
        public static byte Byte(double x)
        {
            return (byte)((x > 255) ? 255 : ((x < 0) ? 0 : x));
        }
        /// <summary>
        /// Преобразовывает случайную величину в тип Byte.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Целое число без знака</returns>
        public static byte Byte(int x)
        {
            return (byte)((x > 255) ? 255 : ((x < 0) ? 0 : x));
        }

        /// <summary>
        /// Преобразовывает случайную величину в тип Byte.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Целое число без знака</returns>
        public static sbyte sByte(double x)
        {
            return (sbyte)((x > 128) ? 128 : ((x < -128) ? -128 : x));
        }
        /// <summary>
        /// Преобразовывает случайную величину в тип Byte.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Целое число без знака</returns>
        public static sbyte sByte(int x)
        {
            return (sbyte)((x > 128) ? 128 : ((x < -128) ? -128 : x));
        }

        /// <summary>
        /// Преобразовывает случайную величину в тип Double.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Double(double x)
        {
            return ((x > 1.0) ? 1.0 : ((x < 0) ? 0 : x));
        }
        /// <summary>
        /// Проверяет лежит ли случайная величина в заданном промежутке.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="xmin">Минимальное значение</param>
        /// <param name="xmax">Максимальное значение</param>
        /// <returns>Логическое значение</returns>
        public static bool IsRange(double x, double xmin, double xmax)
        {
            if (x <= xmax && x >= xmin)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Проверяет лежит ли случайная величина в заданном промежутке.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="xmin">Минимальное значение</param>
        /// <param name="xmax">Максимальное значение</param>
        /// <returns>Логическое значение</returns>
        public static bool IsRange(int x, int xmin, int xmax)
        {
            if (x <= xmax && x >= xmin)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Помещает случайную величину в заданном промежутке.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="xmin">Минимальное значение</param>
        /// <param name="xmax">Максимальное значение</param>
        /// <returns>Логическое значение</returns>
        public static double Range(double x, double xmin, double xmax)
        {
            if (x > xmax)
            {
                return xmax;
            }
            else if (x < xmin)
            {
                return xmin;
            }
            return x;
        }
        /// <summary>
        /// Помещает случайную величину в заданном промежутке.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="xmin">Минимальное значение</param>
        /// <param name="xmax">Максимальное значение</param>
        /// <returns>Логическое значение</returns>
        public static int Range(int x, int xmin, int xmax)
        {
            if (x > xmax)
            {
                return xmax;
            }
            else if (x < xmin)
            {
                return xmin;
            }
            return x;
        }
        #endregion

        #region Singulars
        /// <summary>
        /// Проверяет число на исключение.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Логическое значение</returns>
        public static bool IsSingular(double a)
        {
            if (double.IsNaN(a))
            {
                return true;
            }
            else if (double.IsInfinity(a))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Проверяет комплексное число на исключение.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Логическое значение</returns>
        public static bool IsSingular(Complex a)
        {
            if (IsSingular(a.Re) || IsSingular(a.Im))
            {
                return true;
            }
            return false;
        }
        #endregion

        #region Algebraic
        #region Real number
        /// <summary>
        /// Проверяет является ли число полным квадратом.
        /// </summary>
        /// <param name="n">Целое число</param>
        /// <returns>Логическое значение</returns>
        public static bool IsSquare(double n)
        {
            double sq = (int)Math.Sqrt(n);
            return (sq * sq == n);
        }
        /// <summary>
        /// Проверяет является ли число степенью другого числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="b">Число</param>
        /// <returns>Логическое значение</returns>
        public static bool IsPower(double a, double b)
        {
            double log = Maths.Log(a, b);
            if (IsInteger(log))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Проверяет является ли число целым.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Логическое значение</returns>
        public static bool IsInteger(double a)
        {
            if (a == (int)a)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Проверяет является ли число четным.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Логическое значение</returns>
        public static bool IsEven(double a)
        {
            if (a % 2 == 0)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Проверяет является ли число нечетным.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Логическое значение</returns>
        public static bool IsNotEven(double a)
        {
            return !IsEven(a);
        }
        /// <summary>
        /// Возвращает число, возведенное во вторую степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Pow(double a)
        {
            return Math.Pow(a, 2);
        }
        /// <summary>
        /// Возвращает число, возведенное в степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="b">Степенной показатель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Pow(double a, double b)
        {
            return Math.Pow(a, b);
        }
        /// <summary>
        /// Возвращает экспоненту, возведенную в степень.
        /// </summary>
        /// <param name="a">Степенной показатель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Exp(double a)
        {
            return Math.Pow(E, a);
        }
        /// <summary>
        /// Возвращает натуральный логарифм числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Log(double a)
        {
            return Math.Log(a);
        }
        /// <summary>
        /// Возвращает десятичный логарифм числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Log10(double a)
        {
            return Math.Log(a, 10.0);
        }
        /// <summary>
        /// Возвращает двоичный логарифм числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Log2(double a)
        {
            return Math.Log(a, 2);
        }
        /// <summary>
        /// Возвращает логарифм числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="b">Основание логарифма</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Log(double a, double b)
        {
            return Math.Log(a, b);
        }
        /// <summary>
        /// Возвращает квадратный корень числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Sqrt(double a)
        {
            return Math.Sqrt(a);
        }
        /// <summary>
        /// Возвращает корень числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="b">Степенной показатель корня</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Sqrt(double a, double b)
        {
            return Math.Pow(a, 1.0 / b);
        }
        /// <summary>
        /// Возвращает модуль числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Abs(double a)
        {
            if (a < 0.0)
            {
                return -a;
            }
            return a;
        }
        /// <summary>
        /// Возвращает наибольшее из двух чисел.
        /// </summary>
        /// <param name="a">Первое число</param>
        /// <param name="b">Второе число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Max(double a, double b)
        {
            if (a < b)
            {
                return b;
            }
            return a;
        }
        /// <summary>
        /// Возвращает наибольшее из трех чисел.
        /// </summary>
        /// <param name="a">Первое число</param>
        /// <param name="b">Второе число</param>
        /// <param name="c">Третье число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Max(double a, double b, double c)
        {
            return Max(a, Max(b, c));
        }
        /// <summary>
        /// Возвращает наименьшее из двух чисел.
        /// </summary>
        /// <param name="a">Первое число</param>
        /// <param name="b">Второе число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Min(double a, double b)
        {
            if (a < b)
            {
                return a;
            }
            return b;
        }
        /// <summary>
        /// Возвращает наименьшее из трех чисел.
        /// </summary>
        /// <param name="a">Первое число</param>
        /// <param name="b">Второе число</param>
        /// <param name="c">Третье число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Min(double a, double b, double c)
        {
            return Min(a, Min(b, c));
        }
        /// <summary>
        /// Возвращает знак числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static int Sign(double a)
        {
            if (a < 0)
            {
                return -1;
            }
            else if (a > 0)
            {
                return 1;
            }
            return 0;
        }
        /// <summary>
        /// Возвращает округленное число в меньшую сторону.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Floor(double a)
        {
            return Math.Floor(a);
        }
        /// <summary>
        /// Возвращает округленное число в большую сторону.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Ceil(double a)
        {
            return Math.Ceiling(a);
        }
        /// <summary>
        /// Возвращает округленное число.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Round(double a)
        {
            return Math.Round(a, 0);
        }
        /// <summary>
        /// Возвращает округленное число.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="dig">Количество знаков после запятой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Round(double a, int dig)
        {
            return Math.Round(a, dig);
        }
        /// <summary>
        /// Возвращает число с отброшенной дробной частью.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Fix(double a)
        {
            int ai = (int)a; // Выделение целой части числа
            int sign = Math.Sign(a); // Определение знака
            double c = Maths.Abs(a) - Maths.Abs(ai); // c = |a| - |ai|
            if (c > 0.5)
            {
                return ai + sign; // b = ai + s;
            }
            return ai;
        }
        #endregion

        #region Complex number
        /// <summary>
        /// Возвращает модуль комплексного числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static double Abs(Complex a)
        {
            return a.Abs;
        }
        /// <summary>
        /// Возвращает угол комплексного числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static double Angle(Complex a)
        {
            return a.Angle;
        }
        /// <summary>
        /// Возвращает натуральный логарифм числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Log(Complex a)
        {
            return Log(a, E);
        }
        /// <summary>
        /// Возвращает десятичный логарифм числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Log10(Complex a)
        {
            return Log(a, 10.0);
        }
        /// <summary>
        /// Возвращает двоичный логарифм числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Log2(Complex a)
        {
            return Log(a, 2.0);
        }
        /// <summary>
        /// Возвращает логарифм числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Основание логарифма</param>
        /// <returns>Комплексное число</returns>
        public static Complex Log(Complex a, double b)
        {
            return new Complex(Math.Log(a.Abs), a.Angle) / Math.Log(b);
        }
        /// <summary>
        /// Возвращает экспоненту, возведенную в комплексную степень.
        /// </summary>
        /// <param name="a">Степенной показатель</param>
        /// <returns>Комплексное число</returns>
        public static Complex Exp(Complex a)
        {
            return Pow(E, a);
        }
        /// <summary>
        /// Возвращает число, возведенное в комплексную степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="b">Комплексный показатель степени</param>
        /// <returns>Комплексное число</returns>
        public static Complex Pow(double a, Complex b)
        {
            double r = Math.Pow(a, b.Re);
            return new Complex(r * Math.Cos(b.Im), r * Math.Sin(b.Im));
        }
        /// <summary>
        /// Возвращает число, возведенное в степень.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Показатель степени</param>
        /// <returns>Комплексное число</returns>
        public static Complex Pow(Complex a, double b)
        {
            return Math.Pow(a.Abs, b) * (new Complex(Math.Cos(b * a.Angle), Math.Sin(b * a.Angle)));
        }
        /// <summary>
        /// Возвращает квадратный корень числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Sqrt(Complex a)
        {
            return Maths.Sqrt(a, 2);
        }
        /// <summary>
        /// Возвращает корень числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Показатель степени</param>
        /// <returns>Комплексное число</returns>
        public static Complex Sqrt(Complex a, double b)
        {
            return Maths.FromPolar(Math.Sqrt(a.Abs), a.Angle / b);
        }
        /// <summary>
        /// Возвращает комплексное число.
        /// </summary>
        /// <param name="abs">Модуль комплексного числа</param>
        /// <param name="angle">Угол комплексного числа</param>
        /// <returns>Комплексное число</returns>
        public static Complex FromPolar(double abs, double angle)
        {
            return new Complex(abs * Math.Cos(angle), abs * Math.Sin(angle));
        }
        /// <summary>
        /// Возвращает округленное число.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Round(Complex a)
        {
            return Maths.Round(a, 0);
        }
        /// <summary>
        /// Возвращает округленное число.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="dig">Количество знаков после запятой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Round(Complex a, int dig)
        {
            return new Complex(Math.Round(a.Re, dig), Math.Round(a.Im, dig));
        }
        /// <summary>
        /// Возвращает число с отброшенной дробной частью.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Fix(Complex a)
        {
            return new Complex(Fix(a.Re), Fix(a.Im));
        }
        #endregion
        #endregion

        #region Trigonometric
        #region Real number
        /// <summary>
        /// Возвращает косинус угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Cos(double a)
        {
            return Math.Cos(a);
        }
        /// <summary>
        /// Возвращает синус угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Sin(double a)
        {
            return Math.Sin(a);
        }
        /// <summary>
        /// Возвращает тангенс угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Tg(double a)
        {
            return Math.Sin(a) / Math.Cos(a);
        }
        /// <summary>
        /// Возвращает котангенс угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Ctg(double a)
        {
            return Math.Cos(a) / Math.Sin(a);
        }
        /// <summary>
        /// Возвращает секанс угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Sec(double a)
        {
            return 1.0 / Math.Cos(a);
        }
        /// <summary>
        /// Возвращает косеканс угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Cosc(double a)
        {
            return 1.0 / Math.Sin(a);
        }
        /// <summary>
        /// Возвращает арксинус числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Asin(double a)
        {
            return Math.Asin(a);
        }
        /// <summary>
        /// Возвращает арккосинус числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Acos(double a)
        {
            return Math.Acos(a);
        }
        /// <summary>
        /// Возвращает арктангенс числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Atg(double a)
        {
            return Math.Atan(a);
        }
        /// <summary>
        /// Возвращает арккотангенс числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Actg(double a)
        {
            return (Pi / 2 - Math.Atan(a));
        }
        /// <summary>
        /// Возвращает арксеканс числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Asec(double a)
        {
            return Math.Acos(1.0 / a);
        }
        /// <summary>
        /// Возвращает аркосеканс числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Acosc(double a)
        {
            return Math.Asin(1.0 / a);
        }
        #endregion

        #region Complex number
        /// <summary>
        /// Возвращает косинус угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Cos(Complex a)
        {
            return new Complex(Math.Cos(a.Re) * Math.Cosh(a.Im), -(Math.Sin(a.Re) * Math.Sinh(a.Im)));
        }
        /// <summary>
        /// Возвращает синус угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Sin(Complex a)
        {
            return new Complex(Math.Sin(a.Re) * Math.Cosh(a.Im), Math.Cos(a.Re) * Math.Sinh(a.Im));
        }
        /// <summary>
        /// Возвращает тангенс угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Tg(Complex a)
        {
            return Maths.Sin(a) / Maths.Cos(a);
        }
        /// <summary>
        /// Возвращает котангенс угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Ctg(Complex a)
        {
            return Maths.Cos(a) / Maths.Sin(a);
        }
        /// <summary>
        /// Возвращает секанс угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Sec(Complex a)
        {
            return 1.0 / Maths.Cos(a);
        }
        /// <summary>
        /// Возвращает косеканс угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Cosc(Complex a)
        {
            return 1.0 / Maths.Sin(a);
        }
        /// <summary>
        /// Возвращает арккосинус числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Acos(Complex a)
        {
            return -I * Maths.Log(a + I * Maths.Sqrt(1.0 - a * a));
        }
        /// <summary>
        /// Возвращает арксинус числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Asin(Complex a)
        {
            return -I * Maths.Log(I * a + Maths.Sqrt(1.0 - a * a));
        }
        /// <summary>
        /// Возвращает арктангенс числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Atg(Complex a)
        {
            return I / 2.0 * (Maths.Log(1.0 - I * a) - Maths.Log(1.0 + I * a));
        }
        /// <summary>
        /// Возвращает арккотангенс числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Actg(Complex a)
        {
            return I / 2.0 * (Maths.Log((a - I) / a) - Maths.Log((a + I) / a));
        }
        /// <summary>
        /// Возвращает арксеканс числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Asec(Complex a)
        {
            return Maths.Acos(1.0 / a);
        }
        /// <summary>
        /// Возвращает арккосеканс числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Acosc(Complex a)
        {
            return Maths.Asin(1.0 / a);
        }
        #endregion
        #endregion

        #region Hyperbolic
        #region Real number
        /// <summary>
        /// Возвращает гиперболический синус угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Sh(double a)
        {
            return Math.Sinh(a);
        }
        /// <summary>
        /// Возвращает гиперболический конус угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Ch(double a)
        {
            return Math.Cosh(a);
        }
        /// <summary>
        /// Возвращает гиперболический тангенс угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Th(double a)
        {
            return Math.Sinh(a) / Math.Cosh(a);
        }
        /// <summary>
        /// Возвращает гиперболический котангенс угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Cth(double a)
        {
            return Math.Cosh(a) / Math.Sinh(a);
        }
        /// <summary>
        /// Возвращает гиперболический секанс угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Sch(double a)
        {
            return 1.0 / Math.Cosh(a);
        }
        /// <summary>
        /// Возвращает гиперболический косеканс угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Csch(double a)
        {
            return 1.0 / Math.Sinh(a);
        }
        /// <summary>
        /// Возвращает гиперболический арксинус числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Ash(double a)
        {
            return Math.Log(a + Math.Sqrt(a * a + 1));
        }
        /// <summary>
        /// Возвращает гиперболический арккосинус числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Ach(double a)
        {
            if (a >= 0)
            {
                return Math.Log(a + Math.Sqrt(a * a - 1));
            }
            return 0;
        }
        /// <summary>
        /// Возвращает гиперболический арктангенс числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Ath(double a)
        {
            return 1.0 / 2.0 * Math.Log((1 + a) / (1 - a));
        }
        /// <summary>
        /// Возвращает гиперболический арккотангенс числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Acth(double a)
        {
            return 1.0 / 2.0 * Math.Log((a + 1) / (a - 1));
        }
        /// <summary>
        /// Возвращает гиперболический арксеканс числа.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Asch(double a)
        {
            return Math.Log((1 + Math.Sqrt(1 - a * a)) / a);
        }
        /// <summary>
        /// Возвращает гиперболический арккосеканс числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Acsch(double a)
        {
            if (a < 0)
            {
                return Math.Log((1 - Math.Sqrt(1 + a * a)) / a);
            }
            if (a > 0)
            {
                return Math.Log((1 + Math.Sqrt(1 + a * a)) / a);
            }
            return 0;
        }
        #endregion

        #region Complex number
        /// <summary>
        /// Возвращает гиперболический синус угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Sh(Complex a)
        {
            return new Complex(Math.Sinh(a.Re) * Math.Cos(a.Im), Math.Cosh(a.Re) * Math.Sin(a.Im));
        }
        /// <summary>
        /// Возвращает гиперболический косинус угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Ch(Complex a)
        {
            return new Complex(Math.Cosh(a.Re) * Math.Cos(a.Im), Math.Sinh(a.Re) * Math.Sin(a.Im));
        }
        /// <summary>
        /// Возвращает гиперболический тангенс угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Th(Complex a)
        {
            return Maths.Sh(a) / Maths.Ch(a);
        }
        /// <summary>
        /// Возвращает гиперболический котангенс угла.
        /// </summary>
        /// <param name="a">Угол в радианах</param>
        /// <returns>Комплексное число</returns>
        public static Complex Cth(Complex a)
        {
            return Maths.Ch(a) / Maths.Sh(a);
        }
        /// <summary>
        /// Возвращает гиперболический секанс угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Sch(Complex a)
        {
            return 1.0 / Maths.Ch(a);
        }
        /// <summary>
        /// Возвращает гиперболический косеканс угла.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Csch(Complex a)
        {
            return 1.0 / Maths.Sh(a);
        }
        /// <summary>
        /// Возвращает гиперболический арксинус числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Ash(Complex a)
        {
            return Maths.Log(a + Maths.Sqrt(a * a + 1.0));
        }
        /// <summary>
        /// Возвращает гиперболический арккосинус числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Ach(Complex a)
        {
            return Maths.Log(a + Maths.Sqrt(a * a - 1.0));
        }
        /// <summary>
        /// Возвращает гиперболический арктангенс числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Ath(Complex a)
        {
            return 1.0 / 2.0 * Maths.Log((1.0 + a) / (1.0 - a));
        }
        /// <summary>
        /// Возвращает гиперболический арккотангенс числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Acth(Complex a)
        {
            return 1.0 / 2.0 * Maths.Log((a + 1.0) / (a - 1.0));
        }
        /// <summary>
        /// Возвращает гиперболический арксеканс числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Asch(Complex a)
        {
            return Maths.Log(1.0 / a + Maths.Sqrt(1.0 / a + 1.0) + Maths.Sqrt(1.0 / a - 1.0));
        }
        /// <summary>
        /// Возвращает гиперболический арккосеканс числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex Acsch(Complex a)
        {
            return Maths.Log(1.0 / a + Maths.Sqrt(1.0 / a / a + 1.0));
        }
        #endregion
        #endregion

        #region Modular arithmetic and number theory
        /// <summary>
        /// Проверяет является ли числом простым.
        /// <remarks>
        /// Данный метод основан на переборе всех делителей.
        /// </remarks>
        /// </summary>
        /// <param name="p">Число</param>
        /// <returns>Логическое значение</returns>
        public static bool IsPrime(int p)
        {
            // if number is 2:
            if (p == 2)
            {
                return true;
            }
            // if number is even?
            else if ((p % 2) == 0)
            {
                return false;
            }
            else
            {
                // prime or not?
                int x = Maths.Pollard(p);
                return x == p;
            }
        }
        /// <summary>
        /// Проверяет является ли числом простым.
        /// <remarks>
        /// Данный метод основан на переборе всех делителей.
        /// </remarks>
        /// </summary>
        /// <param name="p">Число</param>
        /// <returns>Логическое значение</returns>
        public static bool IsPrime(long p)
        {
            // if number is 2:
            if (p == 2)
            {
                return true;
            }
            // if number is even?
            else if ((p % 2) == 0)
            {
                return false;
            }
            else
            {
                // prime or not?
                long x = Maths.Pollard(p);
                return x == p;   
            }
        }

        /// <summary>
        /// Возвращает взаимнопростое с "a" число.
        /// </summary>
        /// <param name="a">Целое число</param>
        /// <param name="increment">Начальный вброс</param>
        /// <returns>Целое число со знаком</returns>
        public static int Coprime(int a, int increment = 1)
        {
            int x = 2;
            int p = increment;

            while (x != 1)
            {
                x = Maths.Gcd(a, p);
                p++;
            }

            return p;
        }
        /// <summary>
        /// Возвращает взаимнопростое с "a" число.
        /// </summary>
        /// <param name="a">Целое число</param>
        /// <param name="increment">Начальный вброс</param>
        /// <returns>Целое число со знаком</returns>
        public static long Coprime(long a, long increment = 1)
        {
            long x = 2;
            long p = increment;

            while (x != 1)
            {
                x = Maths.Gcd(a, p);
                p++;
            }

            return p;
        }

        /// <summary>
        /// Возвращает остаток от деления одного числа на другое.
        /// </summary>
        /// <param name="a">Делимое</param>
        /// <param name="n">Делитель</param>
        /// <returns>Целое число со знаком</returns>
        public static int Mod(int a, int n)
        {
            if (n < 0)
                n = -n;

            int r = a % n;
            return r < 0 ? r + n : r;
        }
        /// <summary>
        /// Возвращает остаток от деления одного числа на другое.
        /// </summary>
        /// <param name="a">Делимое</param>
        /// <param name="n">Делитель</param>
        /// <returns>Целое число со знаком</returns>
        public static long Mod(long a, long n)
        {
            if (n < 0)
                n = -n;

            long r = a % n;
            return r < 0 ? r + n : r;
        }
        /// <summary>
        /// Возвращает остаток от деления одного числа на другое.
        /// </summary>
        /// <param name="a">Делимое</param>
        /// <param name="n">Делитель</param>
        /// <returns>Целое число со знаком</returns>
        public static double Mod(double a, double n)
        {
            if (n < 0)
                n = -n;

            double r = a % n;
            return r < 0 ? r + n : r;
        }

        /// <summary>
        /// Возвращает результат возведения числа "a" в степень "x" по модулю p.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="x">Степень</param>
        /// <param name="p">Модуль</param>
        /// <param name="modified">Использовать модифицированный алгоритм или нет</param>
        /// <returns>Целое число со знаком</returns>
        public static int ModPow(int a, int x, int p, bool modified = true)
        {
            if (modified == true)
            {
                return (int)leftmodexp(a, x, p);
            }
            return (int)rightmodexp(a, x, p);
        }
        /// <summary>
        /// Возвращает результат возведения числа "a" в степень "x" по модулю p.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="x">Степень</param>
        /// <param name="p">Модуль</param>
        /// <param name="modified">Использовать модифицированный алгоритм или нет</param>
        /// <returns>Целое число со знаком</returns>
        public static long ModPow(long a, long x, long p, bool modified = true)
        {
            if (modified == true)
            {
                return leftmodexp(a, x, p);
            }
            return rightmodexp(a, x, p);
        }
        /// <summary>
        /// Возвращает результат возведения числа "a" в степень "x" по модулю p.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="x">Степень</param>
        /// <param name="p">Модуль</param>
        /// <returns>Целое число со знаком</returns>
        private static long leftmodexp(long a, long x, long p)
        {
            int[] X = Maths.Decimal2Base(x, 2);
            int t = X.Length, i;
            long y = 1;

            for (i = t - 1; i >= 0; i--)
            {
                y = Maths.Mod(y * y, p);
                if (X[i] == 1)
                {
                    y = Maths.Mod(y * a, p);
                }
            }
            return y;
        }
        /// <summary>
        /// Возвращает результат возведения числа "a" в степень "x" по модулю p.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="x">Степень</param>
        /// <param name="p">Модуль</param>
        /// <returns>Целое число со знаком</returns>
        private static long rightmodexp(long a, long x, long p)
        {
            int[] X = Maths.Decimal2Base(x, 2);
            int t = X.Length, i;
            long y = 1, s = a;

            for (i = 0; i < t; i++)
            {
                if (X[i] == 1)
                {
                    y = Maths.Mod(y * s, p);
                }
                s = Maths.Mod(s * s, p);
            }
            return y;
        }

        /// <summary>
        /// Возвращает обратное число по модулю.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="n">Модуль</param>
        /// <returns>Одномерный массив</returns>
        public static int ModInv(int a, int n)
        {
            int[] U = Euclidean(a, n);
            int gcd = U[0], x = U[1], y = U[2];

            // Условие взаимной простоты:
            if (gcd == 1)
            {
                return (x < 0) ? Maths.Mod(x, n) : x;
            }
            return 0;
        }
        /// <summary>
        /// Возвращает обратное число по модулю.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="n">Модуль</param>
        /// <returns>Одномерный массив</returns>
        public static long ModInv(long a, long n)
        {
            long[] U = Euclidean(a, n);
            long gcd = U[0], x = U[1], y = U[2];

            // Условие взаимной простоты:
            if (gcd == 1)
            {
                return (x < 0) ? Maths.Mod(x, n) : x;
            }
            return 0;
        }

        /// <summary>
        /// Реализует обобщенный алгоритм Евклида.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="n">Модуль</param>
        /// <returns>Одномерный массив</returns>
        public static int[] Euclidean(int a, int n)
        {
            int[] U = new int[3] { a, 1, 0 };
            int[] V = new int[3] { n, 0, 1 };
            int[] T;
            int q;

            while (V[0] != 0)
            {
                q = (int)Maths.Floor(U[0] / V[0]);
                T = new int[3] { Maths.Mod(U[0], V[0]), U[1] - q * V[1], U[2] - q * V[2] };
                U = V;
                V = T;
            }

            return U;
        }
        /// <summary>
        /// Реализует обобщенный алгоритм Евклида.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="n">Модуль</param>
        /// <returns>Одномерный массив</returns>
        public static long[] Euclidean(long a, long n)
        {
            long[] U = new long[3] { a, 1, 0 };
            long[] V = new long[3] { n, 0, 1 };
            long[] T;
            long q;

            while (V[0] != 0)
            {
                q = (long)Maths.Floor(U[0] / V[0]);
                T = new long[3] { Maths.Mod(U[0], V[0]), U[1] - q * V[1], U[2] - q * V[2] };
                U = V;
                V = T;
            }

            return U;
        }

        /// <summary>
        /// Возвращает значение наибольшего общего делителя двух чисел.
        /// </summary>
        /// <param name="a">Целое число</param>
        /// <param name="b">Целое число</param>
        /// <returns>Целое число</returns>
        public static int Gcd(int a, int b)
        {
            int q = Maths.Mod(a, b);
            while (q != 0)
            {
                a = b;
                b = q;
                q = Maths.Mod(a, b);
            }
            return b;
        }
        /// <summary>
        /// Возвращает значение наибольшего общего делителя двух чисел.
        /// </summary>
        /// <param name="a">Целое число</param>
        /// <param name="b">Целое число</param>
        /// <returns>Целое число</returns>
        public static long Gcd(long a, long b)
        {
            long q = Maths.Mod(a, b);
            while (q != 0)
            {
                a = b;
                b = q;
                q = Maths.Mod(a, b);
            }
            return b;
        }

        /// <summary>
        /// Возвращает значение наименьшего общего кратного двух чисел.
        /// </summary>
        /// <param name="a">Целое число</param>
        /// <param name="b">Целое число</param>
        /// <returns>Целое число</returns>
        public static int Lcm(int a, int b)
        {
            return (int)Maths.Abs(a * b) / Gcd(a, b);
        }
        /// <summary>
        /// Возвращает значение наименьшего общего кратного двух чисел.
        /// </summary>
        /// <param name="a">Целое число</param>
        /// <param name="b">Целое число</param>
        /// <returns>Целое число</returns>
        public static long Lcm(long a, long b)
        {
            return (long)Maths.Abs(a * b) / Gcd(a, b);
        }

        /// <summary>
        /// Возвращает массив множителей, из которых состоит число.
        /// </summary>
        /// <param name="n">Целое число</param>
        /// <param name="onlyPrimes">Только простые множители или нет</param>
        /// <returns>Одномерный массив</returns>
        public static int[] Itf(int n, bool onlyPrimes = false)
        {
            int p = n;

            // if collect only prime numbers
            // and "N" includes powers of 2
            if (onlyPrimes)
            {
                int k = 0;

                while (p % 2 == 0)
                {
                    p /= 2;
                    k++;
                }

                if (k > 0)
                    p *= 2;
            }

            // factorization
            var a = new List<int>();
            int div;

            while (p > 1)
            {
                div = Maths.Pollard(p);
                a.Add(div);
                p /= div;
            }

            // distinct or not
            if (onlyPrimes)
            {
                return a.Distinct().ToArray();
            }

            return a.ToArray();
        }
        /// <summary>
        /// Возвращает массив множителей, из которых состоит число.
        /// </summary>
        /// <param name="n">Целое число</param>
        /// <param name="onlyPrimes">Только простые множители или нет</param>
        /// <returns>Одномерный массив</returns>
        public static long[] Itf(long n, bool onlyPrimes = false)
        {
            long p = n;

            // if collect only prime numbers
            // and "N" includes powers of 2
            if (onlyPrimes)
            {
                int k = 0;

                while (p % 2 == 0)
                {
                    p /= 2;
                    k++;
                }

                if (k > 0)
                    p *= 2;
            }

            // factorization
            var a = new List<long>();
            long div;

            while (p > 1)
            {
                div = Maths.Pollard(p);
                a.Add(div);
                p /= div;
            }

            // distinct or not
            if (onlyPrimes)
            {
                return a.Distinct().ToArray();
            }

            return a.ToArray();
        }

        /// <summary>
        /// Возвращает Po-делитель.
        /// </summary>
        /// <param name="n">Целое число</param>
        /// <returns>Целое число</returns>
        public static int Pollard(int n)
        {
            int y = 2, c = 2, x = 2, factor = 1;
            int count;

            while (factor == 1)
            {
                for (count = 1; count <= c && factor <= 1; count++)
                {
                    x = (x * x + 1) % n;
                    factor = Maths.Gcd(x - y, n);
                }

                c *= 2;
                y = x;
            }

            return factor;
        }
        /// <summary>
        /// Возвращает Po-делитель.
        /// </summary>
        /// <param name="n">Целое число</param>
        /// <returns>Целое число</returns>
        public static long Pollard(long n)
        {
            long y = 2, c = 2, x = 2, factor = 1;
            long count;

            while (factor == 1)
            {
                for (count = 1; count <= c && factor <= 1; count++)
                {
                    x = (x * x + 1) % n;
                    factor = Maths.Gcd(x - y, n);
                }

                c *= 2;
                y = x;
            }

            return factor;
        }

        /// <summary>
        /// Возвращает значение функции Эйлера.
        /// </summary>
        /// <param name="n">Натуральное число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static int Etf(int n)
        {
            // factorization with only primes
            int[] itf = Maths.Itf(n, true);
            double radical = 1;
            int length = itf.Length;

            // calculation radical
            for (int i = 0; i < length; i++)
            {
                radical *= 1.0 - 1.0 / itf[i];
            }
            return (int)(n * radical);
        }
        /// <summary>
        /// Возвращает значение функции Эйлера.
        /// </summary>
        /// <param name="n">Натуральное число</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static long Etf(long n)
        {
            // factorization with only primes
            long[] itf = Maths.Itf(n, true);
            double radical = 1;
            int length = itf.Length;

            // calculation radical
            for (int i = 0; i < length; i++)
            {
                radical *= 1.0 - 1.0 / itf[i];
            }
            return (long)(n * radical);
        }

        /// <summary>
        /// Реализует решето для поиска простых чисел.
        /// <remarks>
        /// Рекурсивная реализация оптимизированного по памяти сегментированного решета Эратосфена.
        /// Операционная сложность алгоритма O(N*logN). Сложность по памяти O(Δ), где Δ=sqrt(N).
        /// </remarks>
        /// </summary>
        /// <param name="limit">Число</param>
        /// <returns>Одномерный массив</returns>
        public static int[] Sieve(int limit)
        {
            if (limit <= 2)
            {
                // first prime
                return new int[] { 2 };
            }
            else
            {
                // recursion
                int beta = (int)(Math.Pow(limit, 1.0 / 2)) + 1;
                int[] prime = Sieve(beta);
                bool[] mark;
                int length = prime.Length;
                int start, low, high;
                int i, j, p;
                List<int> list = prime.ToList();

                // do job
                for (low = beta, high = beta + beta; low < limit; low += beta, high += beta)
                {
                    high = Math.Min(high, limit);
                    mark = new bool[beta];

                    for (i = 0; i < length; i++)
                    {
                        p = prime[i];
                        start = (int)((double)low / p) * p;

                        if (start < low)
                            start += p;

                        for (j = start; j < high; j += p)
                        {
                            mark[j - low] = true;
                        }
                    }


                    for (i = low; i < high; i++)
                    {
                        if (!mark[i - low])
                        {
                            list.Add(i);
                        }
                    }
                }

                return list.ToArray();
            }
        }

        /// <summary>
        /// Возвращает радикал целого числа.
        /// </summary>
        /// <param name="n">Число</param>
        /// <returns>Целое число со знаком</returns>
        public static int Radical(int n)
        {
            // factorization
            int[] itf = Maths.Itf(n, true);
            int radical = 1;
            int length = itf.Length;

            // calculation radical
            for (int i = 0; i < length; i++)
            {
                radical *= itf[i];
            }

            return radical;
        }
        /// <summary>
        /// Возвращает радикал целого числа.
        /// </summary>
        /// <param name="n">Число</param>
        /// <returns>Целое число со знаком</returns>
        public static long Radical(long n)
        {
            // factorization
            long[] itf = Maths.Itf(n, true);
            long radical = 1;
            int length = itf.Length;

            // calculation radical
            for (int i = 0; i < length; i++)
            {
                radical *= itf[i];
            }

            return radical;
        }
        #endregion

        #region Private data
        private const int base10 = 10;
        #endregion

        #region Numeral components
        /// <summary>
        /// Возвращает вектор, являющийся представлением десятичного числа в заданной системе счисления.
        /// <remarks>
        /// Пример: 10[10] = {1,0,1,0}[2].
        /// </remarks>
        /// </summary>
        /// <param name="x">Целое число без знака</param>
        /// <param name="newbase">Новое основание числа</param>
        /// <returns>Вектор-строка</returns>
        public static int[] Decimal2Base(long x, int newbase)
        {
            long xc = x;
            int n = NumLength(Math.Abs(xc), newbase);
            int[] X = new int[n];
            int i;

            for (i = 0; i < n; i++)
            {
                X[i] = (int)(Maths.Mod(xc, newbase));
                xc = xc / newbase;
            }

            return X;
        }
        /// <summary>
        /// Возвращает десятичное число, представленное в дисятичной системе счисления.
        /// <remarks>
        /// Пример: {1,0,1,0}[2] = 10[10].
        /// </remarks>
        /// </summary>
        /// <param name="x">Вектор-строка</param>
        /// <param name="thisbase">Основание числа</param>
        /// <returns>Десятичное число</returns>
        public static long Base2Decimal(int[] x, int thisbase)
        {
            int n = x.Length, i;
            long a = 0;

            for (i = 0; i < n; i++)
            {
                a += (long)(x[i] * Maths.Pow(thisbase, i));
            }

            return a;
        }
        /// <summary>
        /// Возвращает число, интерпретирующее заданный вектор в десятичной системе счисления.
        /// <remarks>
        /// Пример: {1,0,1,0}[2] = 1010[10].
        /// </remarks>
        /// </summary>
        /// <param name="x">Вектор-строка</param>
        /// <returns>Десятичное число</returns>
        public static long Vector2Numeral(int[] x)
        {
            int i, n = x.Length;
            long a = 0;

            for (i = 0; i < n; i++)
            {
                a += (long)(x[i] * Maths.Pow(base10, n - i - 1));
            }
            return a;
        }
        /// <summary>
        /// Возвращает вектор-строку, представляющий разбиение десятичного числа на компоненты.
        /// <remarks>
        /// Пример: 1010[10] = {1,0,1,0}[2]
        /// </remarks>
        /// </summary>
        /// <param name="x">Десятичное число</param>
        /// <returns>Вектор-строка</returns>
        public static int[] Numeral2Vector(long x)
        {
            return Decimal2Base(x, base10);
        }
        /// <summary>
        /// Возвращает значение разрядности числа в заданной системе счисления.
        /// </summary>
        /// <param name="x">Целое число без знака</param>
        /// <param name="numbase">Основание</param>
        /// <returns>32-битное число без знаком</returns>
        public static int NumLength(long x, int numbase)
        {
            return (int)Maths.Floor(Maths.Log(x, numbase)) + 1;
        }
        #endregion

        #region Solutions
        /// <summary>
        /// Возвращает значение гипотенузы.
        /// </summary>
        /// <param name="a">Катет</param>
        /// <param name="b">Катет</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Hypotenuse(double a, double b)
        {
            double r = 0.0;
            double absA = System.Math.Abs(a);
            double absB = System.Math.Abs(b);

            if (absA > absB)
            {
                r = b / a;
                r = absA * System.Math.Sqrt(1 + r * r);
            }
            else if (b != 0)
            {
                r = a / b;
                r = absB * System.Math.Sqrt(1 + r * r);
            }

            return r;
        }
        /// <summary>
        /// Реализует решение кубического уравнения вида: 
        /// x^3 + a*x^2 + b*x + c = 0.
        /// </summary>
        /// <param name="a">Коэффициент "a"</param>
        /// <param name="b">Коэффициент "b"</param>
        /// <param name="c">Коэффициент "c"</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Cubic(double a, double b, double c)
        {
            Complex x1 = 0, x2 = 0, x3 = 0;
            double Q = (a * a - 3.0 * b) / 9.0;
            double R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;
            double S = Q * Q * Q - R * R;
            double a3 = a / 3.0;
            double fi, v0, v1;

            // Тригонометрическая формула Виета:
            if (S > 0)
            {
                fi = Math.Acos(R / Math.Pow(Q, 3.0 / 2.0)) / 3.0;
                v0 = -2 * Math.Sqrt(Q);
                v1 = 2.0 / 3 * Math.PI;

                x1 = v0 * Math.Cos(fi) - a3;      //
                x2 = v0 * Math.Cos(fi + v1) - a3; // действительные корни.
                x3 = v0 * Math.Cos(fi - v1) - a3; //
            }
            else if (S < 0)
            {
                if (Q > 0)
                {
                    fi = Maths.Ach(Math.Abs(R) / Math.Pow(Math.Abs(Q), 3.0 / 2.0)) / 3.0;
                    v0 = Math.Sign(R) * Math.Sqrt(Q) * Maths.Ch(fi);
                    v1 = Math.Sqrt(3) * Math.Sqrt(Q) * Maths.Sh(fi);

                    x1 = -2 * v0 - a3;
                    x2 = v0 - a3 + Maths.I * v1;
                    x3 = x2.Conjugate; // комплексно-сопряженный корень.
                }
                else if (Q < 0)
                {
                    fi = Maths.Ash(Math.Abs(R) / Maths.Pow(Math.Abs(Q), 3.0 / 2.0)) / 3.0;
                    v0 = Math.Sign(R) * Math.Sqrt(Math.Abs(Q)) * Maths.Sh(fi);
                    v1 = Math.Sqrt(3) * Math.Sqrt(Math.Abs(Q)) * Maths.Ch(fi); // комплексно-сопряженный корень.

                    x1 = -2 * v0 - a3;
                    x2 = v0 - a3 + Maths.I * v1;
                    x3 = x2.Conjugate; // комплексно-сопряженный корень.
                }
                else if (Q == 0)
                {
                    x1 = -Maths.Sqrt(c - a * a * a / 27.0, 3.0) - a3;
                    v0 = Maths.Abs((a - 3 * x1) * (a + x1) - 4 * b);
                    x2 = Maths.I / 2.0 * Math.Sqrt(v0) - (a + x1) / 2.0;
                    x3 = x2.Conjugate; // комплексно-сопряженный корень.
                }
            }
            else if (S == 0)
            {
                v0 = Math.Pow(R, 1.0 / 3.0);
                x1 = -2 * v0 - a3;
                x2 = x3 = v0 - a3; // пара действительных корней.
            }
            return new Complex[] { x1, x2, x3 };
        }
        /// <summary>
        /// Реализует решение квадратного уравнения вида: 
        /// a*x^2 + b*x + c = 0.
        /// </summary>
        /// <param name="a">Коэффициент "a"</param>
        /// <param name="b">Коэффициент "b"</param>
        /// <param name="c">Коэффициент "c"</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Quadratic(double a, double b, double c)
        {
            // дискриминант уравнения:
            double dis = b * b - 4 * a * c;
            double abs = Math.Sqrt(Math.Abs(dis));
            Complex root = dis < 0 ? new Complex(0, abs) : new Complex(abs, 0);
            Complex q = -0.5 * (b + Math.Sign(b) * root);
            return new Complex[] { q / a, c / q };
        }
        /// <summary>
        /// Реализует решение биквадратного уравнения вида: 
        /// a*x^4 + b*x^2 + c = 0.
        /// </summary>
        /// <param name="a">Коэффициент "a"</param>
        /// <param name="b">Коэффициент "b"</param>
        /// <param name="c">Коэффициент "c"</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] BiQuadratic(double a, double b, double c)
        {
            var s = Quadratic(a, b, c);
            return new Complex[] {     Maths.Sqrt(s[0]), 
                                      -Maths.Sqrt(s[0]), 
                                       Maths.Sqrt(s[1]), 
                                      -Maths.Sqrt(s[1]) };
        }
        #endregion

        #region Givens rotation
        /// <summary>
        /// Реализует построение матрицы вращения Гивенса для пары действительных чисел.
        /// </summary>
        /// <param name="a">Первое число</param>
        /// <param name="b">Второе число</param>
        /// <returns>Матрица</returns>
        public static double[,] Rotation(double a, double b)
        {
            // MATLAB version of
            // Givens rotations:
            double c, s;
            double absx = Maths.Abs(a);

            if (absx == 0)
            {
                c = 0.0;
                s = 1.0;
            }
            else
            {
                double[] v = new double[] { a, b };
                double norm = v.Norm();
                c = absx / norm;
                s = a / absx * (b / norm);
            }

            return new double[,] { { c, s }, { -s, c } };
        }
        /// <summary>
        /// Реализует построение матрицы вращения Гивенса для пары комплексных чисел.
        /// </summary>
        /// <param name="a">Первое число</param>
        /// <param name="b">Второе число</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Rotation(Complex a, Complex b)
        {
            // MATLAB version of
            // Givens rotations:
            Complex c, s;
            Complex absx = Maths.Abs(a);

            if (absx == 0)
            {
                c = 0.0;
                s = 1.0;
            }
            else
            {
                Complex[] v = new Complex[] { a, b };
                double norm = v.Norm();
                c = absx / norm;
                s = a / absx * (b.Conjugate / norm);
            }

            return new Complex[,] { { c, s }, { -s.Conjugate, c } };
        }
        #endregion

        #region Other
        /// <summary>
        /// Нормализует случайную величину относительно {min, max} диапазона.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="min">Минимальное значение, которое принимает случайная величина</param>
        /// <param name="max">Максимальное значение, которое принимает случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static int Normalize(int x, int min, int max)
        {
            int a = max - min;
            int b = x - min;
            int c = (a != 0) ? b / a : x;
            return c;
        }
        /// <summary>
        /// Нормализует случайную величину относительно {min, max} диапазона.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="min">Минимальное значение, которое принимает случайная величина</param>
        /// <param name="max">Максимальное значение, которое принимает случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Normalize(double x, double min, double max)
        {
            double a = max - min;
            double b = x - min;
            double c = (a != 0) ? b / a : x;
            return c;
        }
        #endregion
    }
    /// <summary>
    /// Используется для вычисления расстояний.
    /// </summary>
    public static class Distance
    {
        #region Euclidean distance
        /// <summary>
        /// Возвращает значение расстояния "Euclidean". 
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Euclidean(double[] p, double[] q)
        {
            double sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Math.Pow(p[k] - q[k], 2);
            }

            return Math.Sqrt(sum);
        }
        #endregion

        #region Chebyshev distance
        /// <summary>
        /// Возвращает значение расстояния "Chebyshev". 
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Chebyshev(double[] p, double[] q)
        {
            int n = p.Length;
            double max = Math.Abs(p[0] - q[0]);
            double tmp;

            for (int k = 1; k < n; k++)
            {
                tmp = Math.Abs(p[k] - q[k]);
                max = tmp > max ? tmp : max;
            }

            return max;
        }
        #endregion

        #region Manhattan distance
        /// <summary>
        /// Возвращает значение расстояния "Manhattan". 
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Manhattan(double[] p, double[] q)
        {
            double sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Math.Abs(p[k] - q[k]);
            }

            return sum;
        }
        #endregion

        #region Angular distance
        /// <summary>
        /// Возвращает значение расстояния "Angular". 
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Angular(double[] p, double[] q)
        {
            int n = p.Length;
            double s = 0;
            double x = 0;
            double y = 0;

            for (int i = 0; i < n; i++)
            {
                s += p[i] * q[i];
                x += p[i] * p[i];
                y += q[i] * q[i];
            }

            double den = Math.Sqrt(x) * Math.Sqrt(y);
            double similarity = s == 0 ? 1.0 : 1.0 - (s / den);

            return Math.Acos(similarity);
        }
        #endregion

        #region Bray-Curtis distance
        /// <summary>
        /// Возвращает значение расстояния "Bray-Curtis".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double BrayCurtis(double[] p, double[] q)
        {
            int n = p.Length;
            double x = 0;
            double y = 0;

            for (int i = 0; i < n; i++)
            {
                y += Math.Abs(p[i] - q[i]);
                x += Math.Abs(p[i] + q[i]);
            }

            return y / x;
        }
        #endregion

        #region Canberra distance
        /// <summary>
        /// Возвращает значение расстояния "Canberra".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Canberra(double[] p, double[] q)
        {
            int n = p.Length;
            double sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Math.Abs(p[i] - q[i]) / (Math.Abs(p[i]) + Math.Abs(q[i]));
            }
            return sum;
        }
        #endregion

        #region Dice distance
        /// <summary>
        /// Возвращает значение расстояния "Dice".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Dice(double[] p, double[] q)
        {
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            return (tf + ft) / (double)(2 * tt + ft + tf);
        }
        #endregion

        #region Hellinger distance
        /// <summary>
        /// Возвращает значение расстояния "Hellinger".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Hellinger(double[] p, double[] q)
        {
            int n = p.Length;
            double sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Math.Pow(Math.Sqrt(p[i]) - Math.Sqrt(q[i]), 2);
            }

            return sum / Math.Sqrt(2);
        }
        #endregion

        #region Jaccard distance
        /// <summary>
        /// Возвращает значение расстояния "Jaccard".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Jaccard(double[] p, double[] q)
        {
            int n = p.Length;
            int inter = 0;
            int union = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 || q[i] != 0)
                {
                    if (p[i] == q[i])
                        inter++;
                    union++;
                }
            }

            return (union == 0) ? 0 : 1.0 - (inter / (double)union);
        }
        #endregion

        #region Kulczynski distance
        /// <summary>
        /// Возвращает значение расстояния "Kulczynski".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Kulczynski(double[] p, double[] q)
        {
            // TODO: Rewrite the integer dissimilarities (Yule, Russel-Rao,...)
            // using generics
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            double num = tf + ft - tt + n;
            double den = ft + tf + n;
            return num / den;
        }
        #endregion

        #region Minkowski distance
        /// <summary>
        /// Возвращает значение расстояния "Minkowski".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <param name="order">Порядок</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Minkowski(double[] p, double[] q, double order)
        {
            int n = p.Length;
            double sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Math.Pow(Math.Abs(p[i] - q[i]), order);
            }
            return Math.Pow(sum, 1 / order);
        }
        #endregion

        #region Russel-Rao distance
        /// <summary>
        /// Возвращает значение расстояния "Russel-Rao".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double RusselRao(double[] p, double[] q)
        {
            int n = p.Length;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            return (n - tt) / (double)(n);
        }
        #endregion

        #region Sokal-Michener distance
        /// <summary>
        /// Возвращает значение расстояния "Sokal-Michener".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double SokalMichener(double[] p, double[] q)
        {
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;
            int ff = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] == 1 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] == 1) ft++;
                if (p[i] == 1 && q[i] == 1) tt++;
                if (p[i] == 0 && q[i] == 0) ff++;
            }

            int r = 2 * (tf + ft);
            return r / (double)(ff + tt + r);
        }
        #endregion

        #region Sokal-Sneath distance
        /// <summary>
        /// Возвращает значение расстояния "Sokal-Sneath".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double SokalSneath(double[] p, double[] q)
        {
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            int r = 2 * (tf + ft);
            return r / (double)(tt + r);
        }
        #endregion

        #region Yule distance
        /// <summary>
        /// Возвращает значение расстояния "Yule".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Yule(double[] p, double[] q)
        {
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;
            int ff = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
                if (p[i] == 0 && q[i] == 0) ff++;
            }

            double r = 2 * (tf + ft);
            return r / (tt + ff + r / 2);
        }
        #endregion

        #region Square-Euclidian distance
        /// <summary>
        /// Возвращает значение расстояния "Square-Euclidian".
        /// </summary>
        /// <param name="p">Одномерный массив</param>
        /// <param name="q">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double SquareEuclidian(double[] p, double[] q)
        {
            int n = p.Length;
            double sum = 0.0;
            double u;

            for (int i = 0; i < n; i++)
            {
                u = p[i] - q[i];
                sum += u * u;
            }

            return sum;
        }
        #endregion
    }
    #endregion
}
