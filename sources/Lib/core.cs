// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2020
// Valery Asiryan
// Moscow, Russia

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using System.Xml.Serialization;

namespace UMapx.Core
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                                 UMAPX.CORE
    // **************************************************************************
    // Designed by Valery Asiryan (c), 2015-2020
    // Moscow, Russia.
    // **************************************************************************

    #region Mathematics
    /// <summary>
    /// Uses to implement basic algebraic, trigonometric and hyperbolic operations.
    /// </summary>
    public static class Maths
    {
        #region Constant
        /// <summary>
        /// Exponent.
        /// </summary>
        public const double E = 2.7182818284590452353602874713527;
        /// <summary>
        /// Pi.
        /// </summary>
        public const double Pi = 3.141592653589793238462643383279;
        /// <summary>
        /// Phi (golden number).
        /// </summary>
        public const double Phi = 1.6180339887498948482;
        /// <summary>
        /// Double pi.
        /// </summary>
        public const double Tau = 6.283185307179586476925286766558;
        /// <summary>
        /// Euler-Mascheroni constant.
        /// </summary>
        public const double Gamma = 0.577215664901532860606512090;
        /// <summary>
        /// Square root of number 2.
        /// </summary>
        public const double Sqrt2 = 1.4142135623730950488016887242097;
        /// <summary>
        /// Catalan's constant.
        /// </summary>
        public const double G = 0.915965594177219015054603514932384110774;
        /// <summary>
        /// Apery's constant.
        /// </summary>
        public const double A = 1.202056903159594285399738161511449990764;
        /// <summary>   
        /// Imaginary one.
        /// </summary>
        public static readonly Complex I = Complex.I;
        #endregion

        #region Types and ranges
        /// <summary>
        /// Converts a value to a Byte type.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Byte</returns>
        public static byte Byte(double x)
        {
            return (byte)((x > 255) ? 255 : ((x < 0) ? 0 : x));
        }
        /// <summary>
        /// Converts a value to a Byte type.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Byte</returns>
        public static byte Byte(int x)
        {
            return (byte)((x > 255) ? 255 : ((x < 0) ? 0 : x));
        }

        /// <summary>
        /// Converts a value to a Byte type.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Byte</returns>
        public static sbyte sByte(double x)
        {
            return (sbyte)((x > 128) ? 128 : ((x < -128) ? -128 : x));
        }
        /// <summary>
        /// Converts a value to a Byte type.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Byte</returns>
        public static sbyte sByte(int x)
        {
            return (sbyte)((x > 128) ? 128 : ((x < -128) ? -128 : x));
        }

        /// <summary>
        /// Converts a value to a type Double.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Double precision floating point number</returns>
        public static double Double(double x)
        {
            return ((x > 1.0) ? 1.0 : ((x < 0) ? 0 : x));
        }
        /// <summary>
        /// Checks if value is in the specified range.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="xmax">Maximum value</param>
        /// <returns>Boolean</returns>
        public static bool IsRange(double x, double xmin, double xmax)
        {
            if (x <= xmax && x >= xmin)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks if value is in the specified range.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="xmax">Maximum value</param>
        /// <returns>Boolean</returns>
        public static bool IsRange(int x, int xmin, int xmax)
        {
            if (x <= xmax && x >= xmin)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Crops value in the specified range.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="xmax">Maximum value</param>
        /// <returns>Boolean</returns>
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
        /// Crops value in the specified range.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="xmax">Maximum value</param>
        /// <returns>Boolean</returns>
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
        /// <summary>
        /// Crops value in the specified range.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="xmax">Maximum value</param>
        /// <returns>Boolean</returns>
        public static double Scale(double x, double xmin, double xmax)
        {
            double h = x;

            // bound min
            while (h < xmin)
            {
                h += xmax;
            }

            // bound max
            while (h > xmax)
            {
                h -= xmax;
            }

            return h;
        }
        /// <summary>
        /// Crops value in the specified range.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="xmax">Maximum value</param>
        /// <returns>Boolean</returns>
        public static int Scale(int x, int xmin, int xmax)
        {
            int h = x;

            // bound min
            while (h < xmin)
            {
                h += xmax;
            }

            // bound max
            while (h > xmax)
            {
                h -= xmax;
            }

            return h;
        }
        #endregion

        #region Singulars
        /// <summary>
        /// Checks a number for an exception.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Boolean</returns>
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
        /// Checks a number for an exception.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Boolean</returns>
        public static bool IsSingular(Complex a)
        {
            if (IsSingular(a.Real) || IsSingular(a.Imag))
            {
                return true;
            }
            return false;
        }
        #endregion

        #region Algebraic
        #region Real number
        /// <summary>
        /// Checks if a number is a full square.
        /// </summary>
        /// <param name="n">Integer number</param>
        /// <returns>Boolean</returns>
        public static bool IsSquare(double n)
        {
            double sq = (int)Math.Sqrt(n);
            return (sq * sq == n);
        }
        /// <summary>
        /// Checks whether a number is a power of another number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Number</param>
        /// <returns>Boolean</returns>
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
        /// Checks whether a number is an integer.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Boolean</returns>
        public static bool IsInteger(double a)
        {
            if (a == (int)a)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks whether a number is even.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Boolean</returns>
        public static bool IsEven(double a)
        {
            if (a % 2 == 0)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks whether a number is odd.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Boolean</returns>
        public static bool IsNotEven(double a)
        {
            return !IsEven(a);
        }
        /// <summary>
        /// Returns the number raised to the second power.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Pow(double a)
        {
            return Math.Pow(a, 2);
        }
        /// <summary>
        /// Returns the number raised to the power.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Power</param>
        /// <returns>Double precision floating point number</returns>
        public static double Pow(double a, double b)
        {
            return Math.Pow(a, b);
        }
        /// <summary>
        /// Returns the exponent raised to the power.
        /// </summary>
        /// <param name="a">Power</param>
        /// <returns>Double precision floating point number</returns>
        public static double Exp(double a)
        {
            return Math.Pow(E, a);
        }
        /// <summary>
        /// Returns the natural logarithm of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Log(double a)
        {
            return Math.Log(a);
        }
        /// <summary>
        /// Returns the decimal logarithm of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Log10(double a)
        {
            return Math.Log(a, 10.0);
        }
        /// <summary>
        /// Returns the binary logarithm of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Log2(double a)
        {
            return Math.Log(a, 2);
        }
        /// <summary>
        /// Returns the logarithm of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Base</param>
        /// <returns>Double precision floating point number</returns>
        public static double Log(double a, double b)
        {
            return Math.Log(a, b);
        }
        /// <summary>
        /// Returns the square root of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sqrt(double a)
        {
            return Math.Sqrt(a);
        }
        /// <summary>
        /// Returns the root of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Power</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sqrt(double a, double b)
        {
            return Math.Pow(a, 1.0 / b);
        }
        /// <summary>
        /// Returns the modulus of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Abs(double a)
        {
            if (a < 0.0)
            {
                return -a;
            }
            return a;
        }
        /// <summary>
        /// Returns the largest of two numbers.
        /// </summary>
        /// <param name="a">First number</param>
        /// <param name="b">Second number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Max(double a, double b)
        {
            if (a < b)
            {
                return b;
            }
            return a;
        }
        /// <summary>
        /// Returns the largest of three numbers.
        /// </summary>
        /// <param name="a">First number</param>
        /// <param name="b">Second number</param>
        /// <param name="c">Third number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Max(double a, double b, double c)
        {
            return Max(a, Max(b, c));
        }
        /// <summary>
        /// Returns the smallest of two numbers.
        /// </summary>
        /// <param name="a">First number</param>
        /// <param name="b">Second number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Min(double a, double b)
        {
            if (a < b)
            {
                return a;
            }
            return b;
        }
        /// <summary>
        /// Returns the smallest of three numbers.
        /// </summary>
        /// <param name="a">First number</param>
        /// <param name="b">Second number</param>
        /// <param name="c">Third number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Min(double a, double b, double c)
        {
            return Min(a, Min(b, c));
        }
        /// <summary>
        /// Returns the sign of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns the rounded number down.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Floor(double a)
        {
            return Math.Floor(a);
        }
        /// <summary>
        /// Returns the rounded number up.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Ceil(double a)
        {
            return Math.Ceiling(a);
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Round(double a)
        {
            return Math.Round(a, 0);
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="dig">Digits</param>
        /// <returns>Double precision floating point number</returns>
        public static double Round(double a, int dig)
        {
            return Math.Round(a, dig);
        }
        /// <summary>
        /// Returns number with the fractional part discarded.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Fix(double a)
        {
            int ai = (int)a;
            int sign = Math.Sign(a);
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
        /// Returns the modulus of a complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static double Abs(Complex a)
        {
            return a.Abs;
        }
        /// <summary>
        /// Returns the angle of a complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static double Angle(Complex a)
        {
            return a.Angle;
        }
        /// <summary>
        /// Returns the natural logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Log(Complex a)
        {
            return Log(a, E);
        }
        /// <summary>
        /// Returns the decimal logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Log10(Complex a)
        {
            return Log(a, 10.0);
        }
        /// <summary>
        /// Returns the binary logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Log2(Complex a)
        {
            return Log(a, 2.0);
        }
        /// <summary>
        /// Returns the logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Base</param>
        /// <returns>Complex number</returns>
        public static Complex Log(Complex a, double b)
        {
            return new Complex(Math.Log(a.Abs), a.Angle) / Math.Log(b);
        }
        /// <summary>
        /// Returns the exponent raised to a complex degree.
        /// </summary>
        /// <param name="a">Power</param>
        /// <returns>Complex number</returns>
        public static Complex Exp(Complex a)
        {
            return Pow(E, a);
        }
        /// <summary>
        /// Returns the number raised to a complex power.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Power</param>
        /// <returns>Complex number</returns>
        public static Complex Pow(double a, Complex b)
        {
            double r = Math.Pow(a, b.Real);
            return new Complex(r * Math.Cos(b.Imag), r * Math.Sin(b.Imag));
        }
        /// <summary>
        /// Returns the number raised to the power.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Power</param>
        /// <returns>Complex number</returns>
        public static Complex Pow(Complex a, double b)
        {
            return Math.Pow(a.Abs, b) * (new Complex(Math.Cos(b * a.Angle), Math.Sin(b * a.Angle)));
        }
        /// <summary>
        /// Returns the square root of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Sqrt(Complex a)
        {
            return Maths.Sqrt(a, 2);
        }
        /// <summary>
        /// Returns the root of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Power</param>
        /// <returns>Complex number</returns>
        public static Complex Sqrt(Complex a, double b)
        {
            return Maths.FromPolar(Math.Sqrt(a.Abs), a.Angle / b);
        }
        /// <summary>
        /// Returns complex number.
        /// </summary>
        /// <param name="abs">Module</param>
        /// <param name="angle">Angle</param>
        /// <returns>Complex number</returns>
        public static Complex FromPolar(double abs, double angle)
        {
            return new Complex(abs * Math.Cos(angle), abs * Math.Sin(angle));
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Round(Complex a)
        {
            return Maths.Round(a, 0);
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="dig">Digits</param>
        /// <returns>Complex number</returns>
        public static Complex Round(Complex a, int dig)
        {
            return new Complex(Math.Round(a.Real, dig), Math.Round(a.Imag, dig));
        }
        /// <summary>
        /// Returns number with the fractional part discarded.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Fix(Complex a)
        {
            return new Complex(Fix(a.Real), Fix(a.Imag));
        }
        #endregion
        #endregion

        #region Trigonometric
        #region Real number
        /// <summary>
        /// Returns the cosine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Cos(double a)
        {
            return Math.Cos(a);
        }
        /// <summary>
        /// Returns the sine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sin(double a)
        {
            return Math.Sin(a);
        }
        /// <summary>
        /// Returns the tangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Tg(double a)
        {
            return Math.Sin(a) / Math.Cos(a);
        }
        /// <summary>
        /// Returns the cotangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Ctg(double a)
        {
            return Math.Cos(a) / Math.Sin(a);
        }
        /// <summary>
        /// Returns the secant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sec(double a)
        {
            return 1.0 / Math.Cos(a);
        }
        /// <summary>
        /// Returns the cosecant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Cosc(double a)
        {
            return 1.0 / Math.Sin(a);
        }
        /// <summary>
        /// Returns the arcsine of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Asin(double a)
        {
            return Math.Asin(a);
        }
        /// <summary>
        /// Returns the arccosine of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Acos(double a)
        {
            return Math.Acos(a);
        }
        /// <summary>
        /// Returns the arctangent of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Atg(double a)
        {
            return Math.Atan(a);
        }
        /// <summary>
        /// Returns the arccotangent of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Actg(double a)
        {
            return (Pi / 2 - Math.Atan(a));
        }
        /// <summary>
        /// Returns the arcsecance of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Asec(double a)
        {
            return Math.Acos(1.0 / a);
        }
        /// <summary>
        /// Returns the arccosecant of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Acosc(double a)
        {
            return Math.Asin(1.0 / a);
        }
        #endregion

        #region Complex number
        /// <summary>
        /// Returns the cosine of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Cos(Complex a)
        {
            return new Complex(Math.Cos(a.Real) * Math.Cosh(a.Imag), -(Math.Sin(a.Real) * Math.Sinh(a.Imag)));
        }
        /// <summary>
        /// Returns the sine of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Sin(Complex a)
        {
            return new Complex(Math.Sin(a.Real) * Math.Cosh(a.Imag), Math.Cos(a.Real) * Math.Sinh(a.Imag));
        }
        /// <summary>
        /// Returns the tangent of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Tg(Complex a)
        {
            return Maths.Sin(a) / Maths.Cos(a);
        }
        /// <summary>
        /// Returns the cotangent of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Ctg(Complex a)
        {
            return Maths.Cos(a) / Maths.Sin(a);
        }
        /// <summary>
        /// Returns the secant of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Sec(Complex a)
        {
            return 1.0 / Maths.Cos(a);
        }
        /// <summary>
        /// Returns the cosecant of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Cosc(Complex a)
        {
            return 1.0 / Maths.Sin(a);
        }
        /// <summary>
        /// Returns the arccosine of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Acos(Complex a)
        {
            return -I * Maths.Log(a + I * Maths.Sqrt(1.0 - a * a));
        }
        /// <summary>
        /// Returns the arcsine of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Asin(Complex a)
        {
            return -I * Maths.Log(I * a + Maths.Sqrt(1.0 - a * a));
        }
        /// <summary>
        /// Returns the arctangent of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Atg(Complex a)
        {
            return I / 2.0 * (Maths.Log(1.0 - I * a) - Maths.Log(1.0 + I * a));
        }
        /// <summary>
        /// Returns the arccotangent of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Actg(Complex a)
        {
            return I / 2.0 * (Maths.Log((a - I) / a) - Maths.Log((a + I) / a));
        }
        /// <summary>
        /// Returns the arcsecance of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Asec(Complex a)
        {
            return Maths.Acos(1.0 / a);
        }
        /// <summary>
        /// Returns the arccosecant of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Acosc(Complex a)
        {
            return Maths.Asin(1.0 / a);
        }
        #endregion
        #endregion

        #region Hyperbolic
        #region Real number
        /// <summary>
        /// Returns the hyperbolic sine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sh(double a)
        {
            return Math.Sinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cosine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Ch(double a)
        {
            return Math.Cosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic tangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Th(double a)
        {
            return Math.Sinh(a) / Math.Cosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cotangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Cth(double a)
        {
            return Math.Cosh(a) / Math.Sinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic secant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sch(double a)
        {
            return 1.0 / Math.Cosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cosecant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Csch(double a)
        {
            return 1.0 / Math.Sinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic arcsine of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Ash(double a)
        {
            return Math.Log(a + Math.Sqrt(a * a + 1));
        }
        /// <summary>
        /// Returns the hyperbolic arccosine of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Ach(double a)
        {
            if (a >= 0)
            {
                return Math.Log(a + Math.Sqrt(a * a - 1));
            }
            return 0;
        }
        /// <summary>
        /// Returns the hyperbolic arctangent of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Ath(double a)
        {
            return 1.0 / 2.0 * Math.Log((1 + a) / (1 - a));
        }
        /// <summary>
        /// Returns the hyperbolic arccotangent of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
        public static double Acth(double a)
        {
            return 1.0 / 2.0 * Math.Log((a + 1) / (a - 1));
        }
        /// <summary>
        /// Returns the hyperbolic arcsecance of a number.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Double precision floating point number</returns>
        public static double Asch(double a)
        {
            return Math.Log((1 + Math.Sqrt(1 - a * a)) / a);
        }
        /// <summary>
        /// Returns the hyperbolic arccosecant of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns the hyperbolic sine of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Sh(Complex a)
        {
            return new Complex(Math.Sinh(a.Real) * Math.Cos(a.Imag), Math.Cosh(a.Real) * Math.Sin(a.Imag));
        }
        /// <summary>
        /// Returns the hyperbolic cosine of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Ch(Complex a)
        {
            return new Complex(Math.Cosh(a.Real) * Math.Cos(a.Imag), Math.Sinh(a.Real) * Math.Sin(a.Imag));
        }
        /// <summary>
        /// Returns the hyperbolic tangent of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Th(Complex a)
        {
            return Maths.Sh(a) / Maths.Ch(a);
        }
        /// <summary>
        /// Returns the hyperbolic cotangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Complex number</returns>
        public static Complex Cth(Complex a)
        {
            return Maths.Ch(a) / Maths.Sh(a);
        }
        /// <summary>
        /// Returns the hyperbolic secant of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Sch(Complex a)
        {
            return 1.0 / Maths.Ch(a);
        }
        /// <summary>
        /// Returns the hyperbolic cosecant of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Csch(Complex a)
        {
            return 1.0 / Maths.Sh(a);
        }
        /// <summary>
        /// Returns the hyperbolic arcsine of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Ash(Complex a)
        {
            return Maths.Log(a + Maths.Sqrt(a * a + 1.0));
        }
        /// <summary>
        /// Returns the hyperbolic arccosine of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Ach(Complex a)
        {
            return Maths.Log(a + Maths.Sqrt(a * a - 1.0));
        }
        /// <summary>
        /// Returns the hyperbolic arctangent of a number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <returns>Complex number</returns>
        public static Complex Ath(Complex a)
        {
            return 1.0 / 2.0 * Maths.Log((1.0 + a) / (1.0 - a));
        }
        /// <summary>
        /// Returns the hyperbolic arccotangent of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Acth(Complex a)
        {
            return 1.0 / 2.0 * Maths.Log((a + 1.0) / (a - 1.0));
        }
        /// <summary>
        /// Returns the hyperbolic arcsecance of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Asch(Complex a)
        {
            return Maths.Log(1.0 / a + Maths.Sqrt(1.0 / a + 1.0) + Maths.Sqrt(1.0 / a - 1.0));
        }
        /// <summary>
        /// Returns the hyperbolic arccosecant of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex Acsch(Complex a)
        {
            return Maths.Log(1.0 / a + Maths.Sqrt(1.0 / a / a + 1.0));
        }
        #endregion
        #endregion

        #region Modular arithmetic and number theory
        /// <summary>
        /// Checks if number is prime.
        /// <remarks>
        /// This method is based on enumerating all the divisors.
        /// </remarks>
        /// </summary>
        /// <param name="p">Number</param>
        /// <returns>Boolean</returns>
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
        /// Checks if number is prime.
        /// <remarks>
        /// This method is based on enumerating all the divisors.
        /// </remarks>
        /// </summary>
        /// <param name="p">Number</param>
        /// <returns>Boolean</returns>
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
        /// Returns coprime number.
        /// </summary>
        /// <param name="a">Integer number</param>
        /// <param name="increment">Increment</param>
        /// <returns>Integer number</returns>
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
        /// Returns coprime number.
        /// </summary>
        /// <param name="a">Integer number</param>
        /// <param name="increment">Increment</param>
        /// <returns>Integer number</returns>
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
        /// Returns the remainder of dividing one number by another.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="n">Modulo</param>
        /// <returns>Integer number</returns>
        public static int Mod(int a, int n)
        {
            if (n < 0)
                n = -n;

            int r = a % n;
            return r < 0 ? r + n : r;
        }
        /// <summary>
        /// Returns the remainder of dividing one number by another.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="n">Modulo</param>
        /// <returns>Integer number</returns>
        public static long Mod(long a, long n)
        {
            if (n < 0)
                n = -n;

            long r = a % n;
            return r < 0 ? r + n : r;
        }
        /// <summary>
        /// Returns the remainder of dividing one number by another.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="n">Modulo</param>
        /// <returns>Integer number</returns>
        public static double Mod(double a, double n)
        {
            if (n < 0)
                n = -n;

            double r = a % n;
            return r < 0 ? r + n : r;
        }

        /// <summary>
        /// Returns the result of raising the number "a" to the power of "x" modulo p.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="x">Power</param>
        /// <param name="p">Modulo</param>
        /// <param name="modified">Use modified algorithm or not</param>
        /// <returns>Integer number</returns>
        public static int ModPow(int a, int x, int p, bool modified = true)
        {
            if (modified == true)
            {
                return (int)leftmodexp(a, x, p);
            }
            return (int)rightmodexp(a, x, p);
        }
        /// <summary>
        /// Returns the result of raising the number "a" to the power of "x" modulo p.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="x">Power</param>
        /// <param name="p">Modulo</param>
        /// <param name="modified">Use modified algorithm or not</param>
        /// <returns>Integer number</returns>
        public static long ModPow(long a, long x, long p, bool modified = true)
        {
            if (modified == true)
            {
                return leftmodexp(a, x, p);
            }
            return rightmodexp(a, x, p);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <param name="p"></param>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="x"></param>
        /// <param name="p"></param>
        /// <returns></returns>
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
        /// Returns the inverse number modulo.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="n">Modulo</param>
        /// <returns>Integer number</returns>
        public static int ModInv(int a, int n)
        {
            int[] U = Euclidean(a, n);
            int gcd = U[0], x = U[1], y = U[2];

            if (gcd == 1)
            {
                return (x < 0) ? Maths.Mod(x, n) : x;
            }
            return 0;
        }
        /// <summary>
        /// Returns the inverse number modulo.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="n">Modulo</param>
        /// <returns>Integer number</returns>
        public static long ModInv(long a, long n)
        {
            long[] U = Euclidean(a, n);
            long gcd = U[0], x = U[1], y = U[2];

            if (gcd == 1)
            {
                return (x < 0) ? Maths.Mod(x, n) : x;
            }
            return 0;
        }

        /// <summary>
        /// Implements a generalized Euclidean algorithm.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="n">Modulo</param>
        /// <returns>Array</returns>
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
        /// Implements a generalized Euclidean algorithm.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="n">Modulo</param>
        /// <returns>Array</returns>
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
        /// Returns the value of the greatest common divisor of two numbers.
        /// </summary>
        /// <param name="a">Integer number</param>
        /// <param name="b">Integer number</param>
        /// <returns>Integer number</returns>
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
        /// Returns the value of the greatest common divisor of two numbers.
        /// </summary>
        /// <param name="a">Integer number</param>
        /// <param name="b">Integer number</param>
        /// <returns>Integer number</returns>
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
        /// Returns the value of the least common multiple of two numbers.
        /// </summary>
        /// <param name="a">Integer number</param>
        /// <param name="b">Integer number</param>
        /// <returns>Integer number</returns>
        public static int Lcm(int a, int b)
        {
            return (int)Maths.Abs(a * b) / Gcd(a, b);
        }
        /// <summary>
        /// Returns the value of the least common multiple of two numbers.
        /// </summary>
        /// <param name="a">Integer number</param>
        /// <param name="b">Integer number</param>
        /// <returns>Integer number</returns>
        public static long Lcm(long a, long b)
        {
            return (long)Maths.Abs(a * b) / Gcd(a, b);
        }

        /// <summary>
        /// Returns an array of factors that number consists of.
        /// </summary>
        /// <param name="n">Integer number</param>
        /// <param name="onlyPrimes">Only prime factors or not</param>
        /// <returns>Array</returns>
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
        /// Returns an array of factors that number consists of.
        /// </summary>
        /// <param name="n">Integer number</param>
        /// <param name="onlyPrimes">Only prime factors or not</param>
        /// <returns>Array</returns>
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
        /// Returns the P0-divider.
        /// </summary>
        /// <param name="n">Integer number</param>
        /// <returns>Integer number</returns>
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
        /// Returns the P0-divider.
        /// </summary>
        /// <param name="n">Integer number</param>
        /// <returns>Integer number</returns>
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
        /// Returns the value of the Euler function.
        /// </summary>
        /// <param name="n">Number</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns the value of the Euler function.
        /// </summary>
        /// <param name="n">Number</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Implements a sieve for finding prime numbers.
        /// <remarks>
        /// Recursive implementation of a memory-optimized segmented sieve of Eratosthenes. 
        /// The operational complexity of the O(N* logN) algorithm.The memory complexity is O(Δ), where Δ = sqrt(N).
        /// </remarks>
        /// </summary>
        /// <param name="limit">Number</param>
        /// <returns>Array</returns>
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
        /// Returns the radical of an integer.
        /// </summary>
        /// <param name="n">Number</param>
        /// <returns>Integer number</returns>
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
        /// Returns the radical of an integer.
        /// </summary>
        /// <param name="n">Number</param>
        /// <returns>Integer number</returns>
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
        /// Returns a vector representing the decimal number in the given number system.
        /// <remarks>
        /// Example: 10[10] = {1,0,1,0}[2].
        /// </remarks>
        /// </summary>
        /// <param name="x">Byte</param>
        /// <param name="newbase">Base</param>
        /// <returns>Array</returns>
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
        /// Returns the decimal Number represented in decimal notation.
        /// <remarks>
        /// Example: {1,0,1,0}[2] = 10[10].
        /// </remarks>
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="thisbase">Base</param>
        /// <returns>Integer number</returns>
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
        /// Returns a number that interprets the specified vector in decimal.
        /// <remarks>
        /// Example: {1,0,1,0}[2] = 1010[10].
        /// </remarks>
        /// </summary>
        /// <param name="x">Array</param>
        /// <returns>Integer number</returns>
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
        /// Returns a vector representing the decomposition of a decimal number into components.
        /// <remarks>
        /// Example: 1010[10] = {1,0,1,0}[2]
        /// </remarks>
        /// </summary>
        /// <param name="x">Number</param>
        /// <returns>Array</returns>
        public static int[] Numeral2Vector(long x)
        {
            return Decimal2Base(x, base10);
        }
        /// <summary>
        /// Returns the value of the digit capacity of a number in the given number system.
        /// </summary>
        /// <param name="x">Byte</param>
        /// <param name="numbase">Base</param>
        /// <returns>Integer number</returns>
        public static int NumLength(long x, int numbase)
        {
            return (int)Maths.Floor(Maths.Log(x, numbase)) + 1;
        }
        #endregion

        #region Solutions
        /// <summary>
        /// Returns the value of the hypotenuse.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Number</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Implements the solution of a cubic equation of the form:
        /// x^3 + a*x^2 + b*x + c = 0.
        /// </summary>
        /// <param name="a">Coefficient "a"</param>
        /// <param name="b">Coefficient "b"</param>
        /// <param name="c">Coefficient "c"</param>
        /// <returns>Array</returns>
        public static Complex[] Cubic(double a, double b, double c)
        {
            Complex x1 = 0, x2 = 0, x3 = 0;
            double Q = (a * a - 3.0 * b) / 9.0;
            double R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;
            double S = Q * Q * Q - R * R;
            double a3 = a / 3.0;
            double fi, v0, v1;

            if (S > 0)
            {
                fi = Math.Acos(R / Math.Pow(Q, 3.0 / 2.0)) / 3.0;
                v0 = -2 * Math.Sqrt(Q);
                v1 = 2.0 / 3 * Math.PI;

                x1 = v0 * Math.Cos(fi) - a3;     
                x2 = v0 * Math.Cos(fi + v1) - a3;
                x3 = v0 * Math.Cos(fi - v1) - a3; 
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
                    x3 = x2.Conjugate;
                }
                else if (Q < 0)
                {
                    fi = Maths.Ash(Math.Abs(R) / Maths.Pow(Math.Abs(Q), 3.0 / 2.0)) / 3.0;
                    v0 = Math.Sign(R) * Math.Sqrt(Math.Abs(Q)) * Maths.Sh(fi);
                    v1 = Math.Sqrt(3) * Math.Sqrt(Math.Abs(Q)) * Maths.Ch(fi);

                    x1 = -2 * v0 - a3;
                    x2 = v0 - a3 + Maths.I * v1;
                    x3 = x2.Conjugate;
                }
                else if (Q == 0)
                {
                    x1 = -Maths.Sqrt(c - a * a * a / 27.0, 3.0) - a3;
                    v0 = Maths.Abs((a - 3 * x1) * (a + x1) - 4 * b);
                    x2 = Maths.I / 2.0 * Math.Sqrt(v0) - (a + x1) / 2.0;
                    x3 = x2.Conjugate;
                }
            }
            else if (S == 0)
            {
                v0 = Math.Pow(R, 1.0 / 3.0);
                x1 = -2 * v0 - a3;
                x2 = x3 = v0 - a3;
            }
            return new Complex[] { x1, x2, x3 };
        }
        /// <summary>
        /// Implements a solution to a quadratic equation of the form: 
        /// a*x^2 + b*x + c = 0.
        /// </summary>
        /// <param name="a">Coefficient "a"</param>
        /// <param name="b">Coefficient "b"</param>
        /// <param name="c">Coefficient "c"</param>
        /// <returns>Array</returns>
        public static Complex[] Quadratic(double a, double b, double c)
        {
            double dis = b * b - 4 * a * c;
            double abs = Math.Sqrt(Math.Abs(dis));
            Complex root = dis < 0 ? new Complex(0, abs) : new Complex(abs, 0);
            Complex q = -0.5 * (b + Math.Sign(b) * root);
            return new Complex[] { q / a, c / q };
        }
        /// <summary>
        /// Implements the solution of a biquadratic equation of the form:
        /// a*x^4 + b*x^2 + c = 0.
        /// </summary>
        /// <param name="a">Coefficient "a"</param>
        /// <param name="b">Coefficient "b"</param>
        /// <param name="c">Coefficient "c"</param>
        /// <returns>Array</returns>
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
        /// Implements the construction of the Givens rotation matrix for a pair of real numbers.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Number</param>
        /// <returns>Matrix</returns>
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
        /// Implements the construction of the Givens rotation matrix for a pair of real numbers.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Number</param>
        /// <returns>Matrix</returns>
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
        /// Normalizes a variable relative to the {min, max} range.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="min">Minimum value</param>
        /// <param name="max">Maximum value</param>
        /// <returns>Double precision floating point number</returns>
        public static int Normalize(int x, int min, int max)
        {
            int a = max - min;
            int b = x - min;
            int c = (a != 0) ? b / a : x;
            return c;
        }
        /// <summary>
        /// Normalizes a variable relative to the {min, max} range.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="min">Minimum value</param>
        /// <param name="max">Maximum value</param>
        /// <returns>Double precision floating point number</returns>
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
    /// Uses to calculate distances.
    /// </summary>
    public static class Distance
    {
        #region Euclidean distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value".
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <param name="order">Order</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
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

    #region Serialization
    /// <summary>
    /// Uses for binary serialization of objects.
    /// </summary>
    public static class Binary
    {
        #region Binaries
        /// <summary>
        /// Save data from the file.
        /// </summary>
        /// <param name="stream">Stream</param>
        /// <param name="o">Object</param>
        public static void Save(Stream stream, object o)
        {
            IFormatter bin = new BinaryFormatter();
            bin.Serialize(stream, o);
        }
        /// <summary>
        /// Save data from the file.
        /// </summary>
        /// <param name="fileName">File name</param>
        /// <param name="o">Object</param>
        public static void Save(string fileName, object o)
        {
            FileStream stream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            IFormatter bin = new BinaryFormatter();
            bin.Serialize(stream, o);
            stream.Close();
            stream.Dispose();
        }
        /// <summary>
        /// Load data from the file.
        /// </summary>
        /// <param name="stream">Stream</param>
        public static object Load(Stream stream)
        {
            IFormatter bin = new BinaryFormatter();
            return bin.Deserialize(stream);
        }
        /// <summary>
        /// Load data from the file.
        /// </summary>
        /// <param name="fileName">File name</param>
        public static object Load(string fileName)
        {
            FileStream stream = new FileStream(fileName, FileMode.Open, FileAccess.Read);
            IFormatter bin = new BinaryFormatter();
            object graph = bin.Deserialize(stream);
            stream.Close();
            stream.Dispose();
            return graph;
        }
        #endregion
    }
    /// <summary>
    /// Uses for Xml serialization of objects.
    /// </summary>
    public static class Xml
    {
        #region Binaries
        /// <summary>
        /// Save data from the file.
        /// </summary>
        /// <param name="stream">Stream</param>
        /// <param name="o">Object</param>
        public static void Save(Stream stream, object o)
        {
            XmlSerializer xml = new XmlSerializer(o.GetType());
            xml.Serialize(stream, o);
        }
        /// <summary>
        /// Save data from the file.
        /// </summary>
        /// <param name="fileName">File name</param>
        /// <param name="o">Object</param>
        public static void Save(string fileName, object o)
        {
            FileStream stream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            XmlSerializer xml = new XmlSerializer(o.GetType());
            xml.Serialize(stream, o);
            stream.Close();
            stream.Dispose();
        }
        /// <summary>
        /// Load data from the file.
        /// </summary>
        /// <param name="stream">Stream</param>
        /// <param name="type">Type</param>
        public static object Load(Stream stream, Type type)
        {
            XmlSerializer xml = new XmlSerializer(type);
            return xml.Deserialize(stream);
        }
        /// <summary>
        /// Load data from the file.
        /// </summary>
        /// <param name="fileName">File name</param>
        /// <param name="type">Type</param>
        public static object Load(string fileName, Type type)
        {
            FileStream stream = new FileStream(fileName, FileMode.Open, FileAccess.Read);
            XmlSerializer xml = new XmlSerializer(type);
            object graph = xml.Deserialize(stream);
            stream.Close();
            stream.Dispose();
            return graph;
        }
        #endregion
    }
    #endregion
}
