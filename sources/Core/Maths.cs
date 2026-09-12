using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace UMapx.Core
{
    /// <summary>
    /// Used to implement basic algebraic, trigonometric and hyperbolic operations.
    /// </summary>
    public static partial class Maths
    {
        #region Constant
        /// <summary>
        /// Exponent.
        /// </summary>
        public const float E = 2.7182818284590452353602874713527f;
        /// <summary>
        /// Pi.
        /// </summary>
        public const float Pi = 3.141592653589793238462643383279f;
        /// <summary>
        /// Phi (golden number).
        /// </summary>
        public const float Phi = 1.6180339887498948482f;
        /// <summary>
        /// Double pi.
        /// </summary>
        public const float Tau = 6.283185307179586476925286766558f;
        /// <summary>
        /// Euler-Mascheroni constant.
        /// </summary>
        public const float Gamma = 0.577215664901532860606512090f;
        /// <summary>
        /// Square root of number 2.
        /// </summary>
        public const float Sqrt2 = 1.4142135623730950488016887242097f;
        /// <summary>
        /// Catalan's constant.
        /// </summary>
        public const float G = 0.915965594177219015054603514932384110774f;
        /// <summary>
        /// Apery's constant.
        /// </summary>
        public const float A = 1.202056903159594285399738161511449990764f;
        /// <summary>   
        /// Imaginary one.
        /// </summary>
        public static readonly Complex32 I = Complex32.I;
        #endregion

        #region Types and ranges
        /// <summary>
        /// Converts a value to a Byte type.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Byte.</returns>
        public static byte Byte(float x)
        {
            return (byte)((x > 255) ? 255 : ((x < 0) ? 0 : x));
        }
        /// <summary>
        /// Converts a value to a Byte type.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Byte.</returns>
        public static byte Byte(int x)
        {
            return (byte)((x > 255) ? 255 : ((x < 0) ? 0 : x));
        }

        /// <summary>
        /// Converts a value to an <see cref="sbyte"/> type and clamps it to the range [-128, 127].
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>sbyte.</returns>
        public static sbyte sByte(float x)
        {
            return (sbyte)((x > 127) ? 127 : ((x < -128) ? -128 : x));
        }
        /// <summary>
        /// Converts a value to an <see cref="sbyte"/> type and clamps it to the range [-128, 127].
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>sbyte.</returns>
        public static sbyte sByte(int x)
        {
            return (sbyte)((x > 127) ? 127 : ((x < -128) ? -128 : x));
        }

        /// <summary>
        /// Converts a value to a type float.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <returns>Value.</returns>
        public static float Float(float x)
        {
            return (x > 1.0f) ? 1.0f : ((x < 0) ? 0 : x);
        }
        /// <summary>
        /// Checks if value is in the specified range.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="xmin">Minimum value.</param>
        /// <param name="xmax">Maximum value.</param>
        /// <returns>Boolean.</returns>
        public static bool IsRange(float x, float xmin, float xmax)
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
        /// <param name="x">Value.</param>
        /// <param name="xmin">Minimum value.</param>
        /// <param name="xmax">Maximum value.</param>
        /// <returns>Boolean.</returns>
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
        /// <param name="x">Value.</param>
        /// <param name="xmin">Minimum value.</param>
        /// <param name="xmax">Maximum value.</param>
        /// <returns>float.</returns>
        public static float Range(float x, float xmin, float xmax)
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
        /// <param name="x">Value.</param>
        /// <param name="xmin">Minimum value.</param>
        /// <param name="xmax">Maximum value.</param>
        /// <returns>int.</returns>
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
        /// Wraps a value into the specified range by cyclically adjusting it.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="xmin">Lower bound of the target interval.</param>
        /// <param name="xmax">Period used for cyclic adjustment and upper bound of the interval.</param>
        /// <returns>float.</returns>
        public static float Scale(float x, float xmin, float xmax)
        {
            float h = x;

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
        /// Wraps a value into the specified range by cyclically adjusting it.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="xmin">Lower bound of the target interval.</param>
        /// <param name="xmax">Period used for cyclic adjustment and upper bound of the interval.</param>
        /// <returns>int.</returns>
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
        /// <param name="a">Value.</param>
        /// <returns>Boolean.</returns>
        public static bool IsSingular(float a)
        {
            if (float.IsNaN(a))
            {
                return true;
            }
            else if (float.IsInfinity(a))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks a number for an exception.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Boolean.</returns>
        public static bool IsSingular(Complex32 a)
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
        /// <param name="n">Integer number.</param>
        /// <returns>Boolean.</returns>
        public static bool IsSquare(float n)
        {
            if (n < 0 || float.IsNaN(n) || float.IsInfinity(n)) return false;
            double s = Math.Sqrt(n);
            long r = (long)Math.Floor(s + 0.5);   // nearest integer
                                                  // compare in double to reduce rounding error; relative tolerance:
            double diff = (double)r * r - (double)n;
            return Math.Abs(diff) <= Math.Max(1.0, Math.Abs(n)) * 1e-7;
        }
        /// <summary>
        /// Checks whether a number is a power of another number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Boolean.</returns>
        public static bool IsPower(float a, float b)
        {
            float log = Maths.Log(a, b);
            if (IsInteger(log))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks whether a number is an integer.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Boolean.</returns>
        public static bool IsInteger(float a)
        {
            if (a == Maths.Round(a))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks whether a number is even.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Boolean.</returns>
        public static bool IsEven(float a)
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
        /// <param name="a">Value.</param>
        /// <returns>Boolean.</returns>
        public static bool IsNotEven(float a)
        {
            return !IsEven(a);
        }
        /// <summary>
        /// Returns the number raised to the second power.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Pow(float a)
        {
            return (float)Math.Pow(a, 2);
        }
        /// <summary>
        /// Returns the number raised to the power.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Power.</param>
        /// <returns>Value.</returns>
        public static float Pow(float a, float b)
        {
            return (float)Math.Pow(a, b);
        }
        /// <summary>
        /// Returns the exponent raised to the power.
        /// </summary>
        /// <param name="a">Power.</param>
        /// <returns>Value.</returns>
        public static float Exp(float a)
        {
            return (float)Math.Pow(E, a);
        }
        /// <summary>
        /// Returns the natural logarithm of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Log(float a)
        {
            return (float)Math.Log(a);
        }
        /// <summary>
        /// Returns the decimal logarithm of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Log10(float a)
        {
            return (float)Math.Log(a, 10.0f);
        }
        /// <summary>
        /// Returns the binary logarithm of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Log2(float a)
        {
            return (float)Math.Log(a, 2);
        }
        /// <summary>
        /// Returns the logarithm of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Base.</param>
        /// <returns>Value.</returns>
        public static float Log(float a, float b)
        {
            return (float)Math.Log(a, b);
        }
        /// <summary>
        /// Returns the square root of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Sqrt(float a)
        {
            return (float)Math.Sqrt(a);
        }
        /// <summary>
        /// Returns the root of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Power.</param>
        /// <returns>Value.</returns>
        public static float Sqrt(float a, float b)
        {
            return (float)Math.Pow(a, 1.0f / b);
        }
        /// <summary>
        /// Returns the modulus of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Abs(float a)
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
        /// <param name="a">First number.</param>
        /// <param name="b">Second number.</param>
        /// <returns>Value.</returns>
        public static float Max(float a, float b)
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
        /// <param name="a">First number.</param>
        /// <param name="b">Second number.</param>
        /// <param name="c">Third number.</param>
        /// <returns>Value.</returns>
        public static float Max(float a, float b, float c)
        {
            return Max(a, Max(b, c));
        }
        /// <summary>
        /// Returns the smallest of two numbers.
        /// </summary>
        /// <param name="a">First number.</param>
        /// <param name="b">Second number.</param>
        /// <returns>Value.</returns>
        public static float Min(float a, float b)
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
        /// <param name="a">First number.</param>
        /// <param name="b">Second number.</param>
        /// <param name="c">Third number.</param>
        /// <returns>Value.</returns>
        public static float Min(float a, float b, float c)
        {
            return Min(a, Min(b, c));
        }
        /// <summary>
        /// Returns the sign of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static int Sign(float a)
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
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Floor(float a)
        {
            return (float)Math.Floor(a);
        }
        /// <summary>
        /// Returns the rounded number up.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Ceil(float a)
        {
            return (float)Math.Ceiling(a);
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Round(float a)
        {
            return (float)Math.Round(a, 0);
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="dig">Digits.</param>
        /// <returns>Value.</returns>
        public static float Round(float a, int dig)
        {
            return (float)Math.Round(a, dig);
        }
        #endregion

        #region Complex number
        /// <summary>
        /// Returns the modulus of a complex number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Value.</returns>
        public static float Abs(Complex32 a)
        {
            return a.Abs;
        }
        /// <summary>
        /// Returns the angle of a complex number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Value.</returns>
        public static float Angle(Complex32 a)
        {
            return a.Angle;
        }
        /// <summary>
        /// Returns the natural logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Log(Complex32 a)
        {
            return Complex.Log(a);
        }
        /// <summary>
        /// Returns the decimal logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Log10(Complex32 a)
        {
            return Log(a, 10.0f);
        }
        /// <summary>
        /// Returns the binary logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Log2(Complex32 a)
        {
            return Log(a, 2.0f);
        }
        /// <summary>
        /// Returns the logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <param name="b">Base.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Log(Complex32 a, float b)
        {
            return Maths.Log(a) / Maths.Log(b);
        }
        /// <summary>
        /// Returns the exponent raised to a complex degree.
        /// </summary>
        /// <param name="a">Power.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Exp(Complex32 a)
        {
            float ex = Maths.Exp(a.Real);
            return new Complex32(ex * Maths.Cos(a.Imag),
                                 ex * Maths.Sin(a.Imag));
        }
        /// <summary>
        /// Returns the number raised to a complex power.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <param name="b">Power.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Pow(float a, Complex32 b)
        {
            return Complex.Pow(new Complex(a, 0), b);
        }
        /// <summary>
        /// Returns the number raised to the power.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <param name="b">Power.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Pow(Complex32 a, float b)
        {
            float r = Maths.Pow(a.Abs, b);
            float theta = b * a.Angle;
            return new Complex32(r * Maths.Cos(theta),
                                 r * Maths.Sin(theta));
        }
        /// <summary>
        /// Returns the number raised to the power.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <param name="b">Power.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Pow(Complex32 a, Complex32 b)
        {
            return Maths.Exp(b * Maths.Log(a));
        }
        /// <summary>
        /// Returns the square root of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Sqrt(Complex32 a)
        {
            return Maths.Sqrt(a, 2);
        }
        /// <summary>
        /// Returns the root of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <param name="b">Power.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Sqrt(Complex32 a, float b)
        {
            return Maths.FromPolar(Maths.Pow(a.Abs, 1f / b), a.Angle / b);
        }
        /// <summary>
        /// Returns the root of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <param name="b">Power.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Sqrt(Complex32 a, Complex32 b)
        {
            return Maths.Exp(Maths.Log(a) / b);
        }
        /// <summary>
        /// Returns complex number.
        /// </summary>
        /// <param name="abs">Modulus.</param>
        /// <param name="angle">Angle.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 FromPolar(float abs, float angle)
        {
            return new Complex32(abs * Maths.Cos(angle), abs * Maths.Sin(angle));
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Round(Complex32 a)
        {
            return Maths.Round(a, 0);
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <param name="dig">Digits.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Round(Complex32 a, int dig)
        {
            return new Complex32(Maths.Round(a.Real, dig), Maths.Round(a.Imag, dig));
        }
        #endregion
        #endregion

        #region Trigonometric
        #region Real number
        /// <summary>
        /// Returns the cosine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Cos(float a)
        {
            return (float)Math.Cos(a);
        }
        /// <summary>
        /// Returns the sine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Sin(float a)
        {
            return (float)Math.Sin(a);
        }
        /// <summary>
        /// Returns the tangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Tan(float a)
        {
            return Maths.Sin(a) / Maths.Cos(a);
        }
        /// <summary>
        /// Returns the cotangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Ctan(float a)
        {
            return Maths.Cos(a) / Maths.Sin(a);
        }
        /// <summary>
        /// Returns the secant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Sec(float a)
        {
            return 1.0f / Maths.Cos(a);
        }
        /// <summary>
        /// Returns the cosecant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Cosc(float a)
        {
            return 1.0f / Maths.Sin(a);
        }
        /// <summary>
        /// Returns the arcsine of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Asin(float a)
        {
            return (float)Math.Asin(a);
        }
        /// <summary>
        /// Returns the arccosine of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Acos(float a)
        {
            return (float)Math.Acos(a);
        }
        /// <summary>
        /// Returns the arctangent of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Atan(float a)
        {
            return (float)Math.Atan(a);
        }
        /// <summary>
        /// Returns the arctangent2 of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Value.</returns>
        public static float Atan2(float a, float b)
        {
            return (float)Math.Atan2(a, b);
        }
        /// <summary>
        /// Returns the arccotangent of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Actan(float a)
        {
            return (float)Math.Atan2(1.0, a);
        }
        /// <summary>
        /// Returns the arcsecant of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Asec(float a)
        {
            return (float)Math.Acos(1.0 / a);
        }
        /// <summary>
        /// Returns the arccosecant of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Acosc(float a)
        {
            return (float)Math.Asin(1.0 / a);
        }
        #endregion

        #region Complex number
        /// <summary>
        /// Returns the cosine of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Cos(Complex32 a)
        {
            return new Complex32(Maths.Cos(a.Real) * Maths.Cosh(a.Imag), -(Maths.Sin(a.Real) * Maths.Sinh(a.Imag)));
        }
        /// <summary>
        /// Returns the sine of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Sin(Complex32 a)
        {
            return new Complex32(Maths.Sin(a.Real) * Maths.Cosh(a.Imag), Maths.Cos(a.Real) * Maths.Sinh(a.Imag));
        }
        /// <summary>
        /// Returns the tangent of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Tan(Complex32 a)
        {
            return Maths.Sin(a) / Maths.Cos(a);
        }
        /// <summary>
        /// Returns the cotangent of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Ctan(Complex32 a)
        {
            return Maths.Cos(a) / Maths.Sin(a);
        }
        /// <summary>
        /// Returns the secant of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Sec(Complex32 a)
        {
            return 1.0 / Maths.Cos(a);
        }
        /// <summary>
        /// Returns the cosecant of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Cosc(Complex32 a)
        {
            return 1.0 / Maths.Sin(a);
        }
        /// <summary>
        /// Returns the arccosine of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Acos(Complex32 a)
        {
            return -I * Maths.Log(a + I * Maths.Sqrt(1.0 - a * a));
        }
        /// <summary>
        /// Returns the arcsine of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Asin(Complex32 a)
        {
            return -I * Maths.Log(I * a + Maths.Sqrt(1.0 - a * a));
        }
        /// <summary>
        /// Returns the arctangent of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Atan(Complex32 a)
        {
            return I / 2.0 * (Maths.Log(1.0 - I * a) - Maths.Log(1.0 + I * a));
        }
        /// <summary>
        /// Returns the arctangent2 of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Atan2(Complex32 a, Complex32 b)
        {
            float re = b.Real - a.Imag;
            float im = b.Imag + a.Real;

            float angle = Maths.Atan2(im, re);
            return new Complex32(angle, 0f);
        }
        /// <summary>
        /// Returns the arccotangent of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <remarks>Uses principal atan(1/a), with value pi/2 at zero. On the imaginary cuts,
        /// the real part has the sign of Im(1/a). The real overload uses the interval (0, pi).</remarks>
        /// <returns>Complex number.</returns>
        public static Complex32 Actan(Complex32 a)
        {
            // Principal atan(1/z); use the continuous real-axis value at zero.
            if (a.Real == 0 && a.Imag == 0) return new Complex32((float)(Math.PI / 2), 0);
            return PrincipalAtan(Complex.One / (Complex)a);
        }
        /// <summary>
        /// Returns the arcsecant of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Asec(Complex32 a)
        {
            return Maths.Acos(1.0 / a);
        }
        /// <summary>
        /// Returns the arccosecant of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Acosc(Complex32 a)
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
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Sinh(float a)
        {
            return (float)Math.Sinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cosine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Cosh(float a)
        {
            return (float)Math.Cosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic tangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Tanh(float a)
        {
            return (float)Math.Tanh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cotangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Ctanh(float a)
        {
            return (float)(1.0 / Math.Tanh(a));
        }
        /// <summary>
        /// Returns the hyperbolic secant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Sech(float a)
        {
            return (float)(1.0 / Math.Cosh(a));
        }
        /// <summary>
        /// Returns the hyperbolic cosecant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Cosch(float a)
        {
            return (float)(1.0 / Math.Sinh(a));
        }
        /// <summary>
        /// Returns the hyperbolic arcsine of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Asinh(float a)
        {
            return (float)RealAsinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic arccosine of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Acosh(float a)
        {
            return (float)RealAcosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic arctangent of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Atanh(float a)
        {
            return (float)RealAtanh(a);
        }
        /// <summary>
        /// Returns the hyperbolic arccotangent of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Actanh(float a)
        {
            return (float)RealAtanh(1.0 / a);
        }
        /// <summary>
        /// Returns the hyperbolic arcsecant of a number.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Value.</returns>
        public static float Asech(float a)
        {
            return (float)RealAcosh(1.0 / a);
        }
        /// <summary>
        /// Returns the hyperbolic arccosecant of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Value.</returns>
        public static float Acosch(float a)
        {
            return (float)RealAsinh(1.0 / a);
        }
        #endregion

        #region Complex number
        /// <summary>
        /// Returns the hyperbolic sine of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Sinh(Complex32 a)
        {
            return new Complex32(Maths.Sinh(a.Real) * Maths.Cos(a.Imag), Maths.Cosh(a.Real) * Maths.Sin(a.Imag));
        }
        /// <summary>
        /// Returns the hyperbolic cosine of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Cosh(Complex32 a)
        {
            return new Complex32(Maths.Cosh(a.Real) * Maths.Cos(a.Imag), Maths.Sinh(a.Real) * Maths.Sin(a.Imag));
        }
        /// <summary>
        /// Returns the hyperbolic tangent of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Tanh(Complex32 a)
        {
            if (Math.Abs(a.Real) > 100 && !float.IsNaN(a.Imag) && !float.IsInfinity(a.Imag))
                return new Complex32(a.Real < 0 ? -1 : 1, 0);
            return ComplexSinh(a) / ComplexCosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cotangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Ctanh(Complex32 a)
        {
            if (Math.Abs(a.Real) > 100 && !float.IsNaN(a.Imag) && !float.IsInfinity(a.Imag))
                return new Complex32(a.Real < 0 ? -1 : 1, 0);
            return ComplexCosh(a) / ComplexSinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic secant of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Sech(Complex32 a)
        {
            if (Math.Abs(a.Real) > 105 && !float.IsNaN(a.Imag) && !float.IsInfinity(a.Imag)) return new Complex32(0, 0);
            return Complex.One / ComplexCosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cosecant of an angle.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Cosch(Complex32 a)
        {
            if (Math.Abs(a.Real) > 105 && !float.IsNaN(a.Imag) && !float.IsInfinity(a.Imag)) return new Complex32(0, 0);
            return Complex.One / ComplexSinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic arcsine of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Asinh(Complex32 a)
        {
            if (a.Real < 0 || (a.Real == 0 && a.Imag < 0)) return -Asinh(-a);
            Complex z = a;
            if (z.Magnitude < 1e-4) return a;
            return Complex.Log(z + Complex.Sqrt(z * z + 1));
        }
        /// <summary>
        /// Returns the hyperbolic arccosine of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Acosh(Complex32 a)
        {
            Complex z = a;
            // The product of principal square roots selects the right branch in both half-planes.
            if ((z - 1).Magnitude < 0.5)
                return ComplexLogOnePlus(z - 1 + Complex.Sqrt(z - 1) * Complex.Sqrt(z + 1));
            return Complex.Log(z + Complex.Sqrt(z - 1) * Complex.Sqrt(z + 1));
        }
        /// <summary>
        /// Returns the hyperbolic arctangent of a number.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Atanh(Complex32 a)
        {
            return 1.0 / 2.0 * Maths.Log((1.0 + a) / (1.0 - a));
        }
        /// <summary>
        /// Returns the hyperbolic arccotangent of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Actanh(Complex32 a)
        {
            return 1.0 / 2.0 * Maths.Log((a + 1.0) / (a - 1.0));
        }
        /// <summary>
        /// Returns the hyperbolic arcsecant of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Asech(Complex32 a)
        {
            var inv = 1.0f / a;
            return Maths.Log(inv + Maths.Sqrt(inv + 1.0f) * Maths.Sqrt(inv - 1.0f));
        }
        /// <summary>
        /// Returns the hyperbolic arccosecant of a number.
        /// </summary>
        /// <param name="a">Complex number.</param>
        /// <returns>Complex number.</returns>
        public static Complex32 Acosch(Complex32 a)
        {
            var inv = 1.0f / a;
            return Maths.Log(inv + Maths.Sqrt(inv * inv + 1.0f));
        }
        #endregion
        #endregion

        #region Modular arithmetic and number theory
        /// <summary>
        /// Checks if number is prime.
        /// </summary>
        /// <remarks>
        /// Uses deterministic Miller-Rabin tests over the full signed integer range.
        /// </remarks>
        /// <param name="p">Value.</param>
        /// <returns>Boolean.</returns>
        public static bool IsPrime(int p)
        {
            return p >= 2 && IsPrimeUnsigned((ulong)p);
        }
        /// <summary>
        /// Checks if number is prime.
        /// </summary>
        /// <remarks>
        /// Uses deterministic Miller-Rabin tests over the full signed integer range.
        /// </remarks>
        /// <param name="p">Value.</param>
        /// <returns>Boolean.</returns>
        public static bool IsPrime(long p)
        {
            return p >= 2 && IsPrimeUnsigned((ulong)p);
        }

        /// <summary>
        /// Returns coprime number.
        /// </summary>
        /// <param name="a">Integer number.</param>
        /// <param name="increment">Inclusive starting value for the coprime search.</param>
        /// <returns>Integer number.</returns>
        public static int Coprime(int a, int increment = 1)
        {
            if (a == 0)
            {
                if (increment <= -1) return -1;
                if (increment <= 1) return 1;
                throw new OverflowException("No coprime exists at or above the starting value.");
            }
            while (UnsignedGcd(UnsignedMagnitude(a), UnsignedMagnitude(increment)) != 1)
                increment = checked(increment + 1);
            return increment;
        }
        /// <summary>
        /// Returns coprime number.
        /// </summary>
        /// <param name="a">Integer number.</param>
        /// <param name="increment">Inclusive starting value for the coprime search.</param>
        /// <returns>Integer number.</returns>
        public static long Coprime(long a, long increment = 1)
        {
            if (a == 0)
            {
                if (increment <= -1) return -1;
                if (increment <= 1) return 1;
                throw new OverflowException("No coprime exists at or above the starting value.");
            }
            while (UnsignedGcd(UnsignedMagnitude(a), UnsignedMagnitude(increment)) != 1)
                increment = checked(increment + 1);
            return increment;
        }

        /// <summary>
        /// Returns the remainder of dividing one number by another.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="n">Modulo.</param>
        /// <returns>Integer number.</returns>
        public static int Mod(int a, int n)
        {
            return (int)Mod((long)a, n);
        }
        /// <summary>
        /// Returns the remainder of dividing one number by another.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="n">Modulo.</param>
        /// <returns>Integer number.</returns>
        public static long Mod(long a, long n)
        {
            if (n == -1) return 0; // Includes long.MinValue without signed division overflow.
            long remainder = a % n;
            return remainder < 0 ? (n > 0 ? remainder + n : remainder - n) : remainder;
        }
        /// <summary>
        /// Returns the remainder of dividing one number by another.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="n">Modulo.</param>
        /// <returns>float.</returns>
        public static float Mod(float a, float n)
        {
            if (n < 0)
                n = -n;

            float r = a % n;
            return r < 0 ? r + n : r;
        }

        /// <summary>
        /// Returns the result of raising the number "a" to the power of "x" modulo p.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="x">Power.</param>
        /// <param name="p">Modulo.</param>
        /// <param name="modified">Use modified algorithm or not.</param>
        /// <returns>Integer number.</returns>
        public static int ModPow(int a, int x, int p, bool modified = true)
        {
            return (int)ModPow((long)a, x, p, modified);
        }
        /// <summary>
        /// Returns the result of raising the number "a" to the power of "x" modulo p.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="x">Power.</param>
        /// <param name="p">Modulo.</param>
        /// <param name="modified">Use modified algorithm or not.</param>
        /// <returns>Integer number.</returns>
        public static long ModPow(long a, long x, long p, bool modified = true)
        {
            if (x < 0) throw new ArgumentOutOfRangeException(nameof(x), "The exponent must be nonnegative.");
            if (p == 0) throw new DivideByZeroException();
            return modified ? Leftmodexp(a, x, p) : Rightmodexp(a, x, p);
        }
        /// <summary>
        /// Computes modular exponentiation using the left-to-right binary method.
        /// </summary>
        /// <param name="a">Base value.</param>
        /// <param name="x">Exponent.</param>
        /// <param name="p">Modulus.</param>
        /// <returns>Result of a^x mod p.</returns>
        private static long Leftmodexp(long a, long x, long p)
        {
            ulong modulus = UnsignedMagnitude(p);
            ulong value = (ulong)Mod(a, p), result = 1 % modulus;
            for (int bit = 62; bit >= 0; bit--)
            {
                result = MultiplyModulo(result, result, modulus);
                if (((x >> bit) & 1) != 0) result = MultiplyModulo(result, value, modulus);
            }
            return (long)result;
        }
        /// <summary>
        /// Computes modular exponentiation using the right-to-left binary method.
        /// </summary>
        /// <param name="a">Base value.</param>
        /// <param name="x">Exponent.</param>
        /// <param name="p">Modulus.</param>
        /// <returns>Result of a^x mod p.</returns>
        private static long Rightmodexp(long a, long x, long p)
        {
            ulong modulus = UnsignedMagnitude(p);
            ulong value = (ulong)Mod(a, p), result = 1 % modulus;
            while (x != 0)
            {
                if ((x & 1) != 0) result = MultiplyModulo(result, value, modulus);
                x >>= 1;
                if (x != 0) value = MultiplyModulo(value, value, modulus);
            }
            return (long)result;
        }

        /// <summary>
        /// Returns the inverse number modulo.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="n">Modulo.</param>
        /// <returns>Integer number.</returns>
        public static int ModInv(int a, int n)
        {
            if (n == 0) throw new DivideByZeroException();
            BigInteger[] result = ExtendedGcd(a, n);
            if (result[0] != 1) return 0;
            BigInteger modulus = BigInteger.Abs(n);
            return (int)((result[1] % modulus + modulus) % modulus);
        }
        /// <summary>
        /// Returns the inverse number modulo.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="n">Modulo.</param>
        /// <returns>Integer number.</returns>
        public static long ModInv(long a, long n)
        {
            if (n == 0) throw new DivideByZeroException();
            BigInteger[] result = ExtendedGcd(a, n);
            if (result[0] != 1) return 0;
            BigInteger modulus = BigInteger.Abs(n);
            return (long)((result[1] % modulus + modulus) % modulus);
        }

        /// <summary>
        /// Implements a generalized Euclidean algorithm.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="n">Modulo.</param>
        /// <returns>Array.</returns>
        public static int[] Euclidean(int a, int n)
        {
            BigInteger[] result = ExtendedGcd(a, n);
            return new[] { (int)result[0], (int)result[1], (int)result[2] };
        }
        /// <summary>
        /// Implements a generalized Euclidean algorithm.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="n">Modulo.</param>
        /// <returns>Array.</returns>
        public static long[] Euclidean(long a, long n)
        {
            BigInteger[] result = ExtendedGcd(a, n);
            return new[] { (long)result[0], (long)result[1], (long)result[2] };
        }

        /// <summary>
        /// Returns the value of the greatest common divisor of two numbers.
        /// </summary>
        /// <param name="a">Integer number.</param>
        /// <param name="b">Integer number.</param>
        /// <exception cref="OverflowException">The nonnegative GCD does not fit in Int32.</exception>
        /// <returns>Integer number.</returns>
        public static int Gcd(int a, int b)
        {
            return checked((int)UnsignedGcd(UnsignedMagnitude(a), UnsignedMagnitude(b)));
        }
        /// <summary>
        /// Returns the value of the greatest common divisor of two numbers.
        /// </summary>
        /// <param name="a">Integer number.</param>
        /// <param name="b">Integer number.</param>
        /// <exception cref="OverflowException">The nonnegative GCD does not fit in Int64.</exception>
        /// <returns>Integer number.</returns>
        public static long Gcd(long a, long b)
        {
            return checked((long)UnsignedGcd(UnsignedMagnitude(a), UnsignedMagnitude(b)));
        }

        /// <summary>
        /// Returns the value of the least common multiple of two numbers.
        /// </summary>
        /// <param name="a">Integer number.</param>
        /// <param name="b">Integer number.</param>
        /// <remarks>Returns zero if either input is zero.</remarks>
        /// <exception cref="OverflowException">The nonnegative LCM does not fit in Int32.</exception>
        /// <returns>Integer number.</returns>
        public static int Lcm(int a, int b)
        {
            if (a == 0 || b == 0) return 0;
            BigInteger gcd = BigInteger.GreatestCommonDivisor(a, b);
            return (int)BigInteger.Abs((BigInteger)a / gcd * b);
        }
        /// <summary>
        /// Returns the value of the least common multiple of two numbers.
        /// </summary>
        /// <param name="a">Integer number.</param>
        /// <param name="b">Integer number.</param>
        /// <remarks>Returns zero if either input is zero.</remarks>
        /// <exception cref="OverflowException">The nonnegative LCM does not fit in Int64.</exception>
        /// <returns>Integer number.</returns>
        public static long Lcm(long a, long b)
        {
            if (a == 0 || b == 0) return 0;
            BigInteger gcd = BigInteger.GreatestCommonDivisor(a, b);
            return (long)BigInteger.Abs((BigInteger)a / gcd * b);
        }

        /// <summary>
        /// Returns an array of factors that number consists of.
        /// </summary>
        /// <param name="n">Integer number.</param>
        /// <param name="onlyPrimes">Return distinct prime factors when true; include multiplicities otherwise.</param>
        /// <returns>Array.</returns>
        public static int[] Itf(int n, bool onlyPrimes = false)
        {
            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), "Factorization requires a positive integer.");
            var factors = new List<ulong>();
            FactorInteger((ulong)n, factors);
            factors.Sort();
            return (onlyPrimes ? factors.Distinct() : factors).Select(x => (int)x).ToArray();
        }
        /// <summary>
        /// Returns an array of factors that number consists of.
        /// </summary>
        /// <param name="n">Integer number.</param>
        /// <param name="onlyPrimes">Return distinct prime factors when true; include multiplicities otherwise.</param>
        /// <returns>Array.</returns>
        public static long[] Itf(long n, bool onlyPrimes = false)
        {
            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), "Factorization requires a positive integer.");
            var factors = new List<ulong>();
            FactorInteger((ulong)n, factors);
            factors.Sort();
            return (onlyPrimes ? factors.Distinct() : factors).Select(x => (long)x).ToArray();
        }

        /// <summary>
        /// Returns a proper divisor of a composite positive integer, or the input for a prime or one.
        /// </summary>
        /// <param name="n">Integer number.</param>
        /// <returns>Integer number.</returns>
        public static int Pollard(int n)
        {
            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), "Factorization requires a positive integer.");
            if (n == 1 || IsPrime(n)) return n;
            return (int)FindDivisor((ulong)n);
        }
        /// <summary>
        /// Returns a proper divisor of a composite positive integer, or the input for a prime or one.
        /// </summary>
        /// <param name="n">Integer number.</param>
        /// <returns>Integer number.</returns>
        public static long Pollard(long n)
        {
            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), "Factorization requires a positive integer.");
            if (n == 1 || IsPrime(n)) return n;
            return (long)FindDivisor((ulong)n);
        }

        /// <summary>
        /// Returns the value of the Euler function.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <returns>Value.</returns>
        public static int Etf(int n)
        {
            int result = n;
            foreach (int prime in Itf(n, true)) result = result / prime * (prime - 1);
            return result;
        }
        /// <summary>
        /// Returns the value of the Euler function.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <returns>Value.</returns>
        public static long Etf(long n)
        {
            long result = n;
            foreach (long prime in Itf(n, true)) result = result / prime * (prime - 1);
            return result;
        }

        /// <summary>
        /// Implements a sieve for finding prime numbers.
        /// </summary>
        /// <remarks>
        /// Recursive implementation of a memory-optimized segmented sieve of Eratosthenes. 
        /// The operational complexity of the O(N* logN) algorithm.The memory complexity is O(Δ), where Δ = sqrt(N).
        /// </remarks>
        /// <param name="limit">Value.</param>
        /// <returns>Array.</returns>
        public static int[] Sieve(int limit)
        {
            if (limit < 2) return Array.Empty<int>();
            if (limit == 2) return new[] { 2 };

            // 1) Base primes up to floor(sqrt(limit)) using odd-only sieve
            int sqrt = (int)Math.Sqrt(limit);
            var basePrimes = BuildBasePrimes(sqrt); // includes 2

            // Reserve output capacity using π(n) ~ n/(ln n - 1.08366)
            int cap = limit >= 17 ? (int)(limit / (Math.Log(limit) - 1.08366)) : 8;
            var primes = new List<int>(Math.Max(8, cap))
            {
                2 // we sieve only odds below
            };

            // 2) Segmented sieve over [3..limit], odd-only with bitset
            const int defaultSegmentOddCount = 1 << 20; // how many odds per segment (~128 KiB bitset)
            int segmentOddCount = defaultSegmentOddCount;

            // Pre-allocate bitset for the largest segment: 1 bit per odd
            ulong[] bits = new ulong[(segmentOddCount + 63) >> 6];

            // Use long for exclusive upper bound to avoid overflow at int.MaxValue
            long limitPlus1 = (long)limit + 1;

            // Iterate segments, [low, high), low/high are integers; we mark only odds inside
            for (long low = 3; low <= limit;)
            {
                long high = Math.Min(low + ((long)segmentOddCount << 1), limitPlus1);

                // First odd ≥ low
                long firstOdd = (low | 1);

                // Number of odd integers in [firstOdd, high):
                // count = ceil((high - firstOdd)/2) = (high - firstOdd + 1) >> 1, clamped to ≥0
                int oddCount = high > firstOdd ? (int)((high - firstOdd + 1) >> 1) : 0;

                // Clear only the slice we use
                int words = (oddCount + 63) >> 6;
                if (words > 0) Array.Clear(bits, 0, words);

                // Mark composites using base primes (skip 2; we only store odds)
                for (int t = 0; t < basePrimes.Length; t++)
                {
                    int p = basePrimes[t];
                    if (p == 2) continue;

                    long pp = 1L * p * p;
                    if (pp >= high) break; // nothing to mark in this segment

                    // First multiple of p inside [firstOdd, high)
                    long start = pp > firstOdd ? pp : (firstOdd + p - 1) / p * p;
                    if ((start & 1L) == 0) start += p;        // ensure odd composite
                    int step = p << 1;                         // jump between odd multiples

                    for (long j = start; j < high; j += step)
                    {
                        int idx = (int)((j - firstOdd) >> 1);  // 0 ≤ idx < oddCount
                        bits[idx >> 6] |= 1UL << (idx & 63);
                    }
                }

                // Emit primes from this segment
                for (int i = 0; i < oddCount; i++)
                {
                    if ((bits[i >> 6] & (1UL << (i & 63))) == 0)
                    {
                        int n = (int)(firstOdd + (i << 1));
                        if (n <= limit) primes.Add(n);
                    }
                }

                // Next segment starts at high (exclusive)
                low = high;
            }

            return primes.ToArray();
        }
        /// <summary>
        /// Odd-only sieve up to n (inclusive). Returns base primes incl. 2.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <returns>Array.</returns>
        private static int[] BuildBasePrimes(int n)
        {
            if (n < 2) return Array.Empty<int>();
            if (n == 2) return new[] { 2 };

            // Odd candidates: 3,5,7,...,n
            int m = (n - 1) >> 1; // index i ↦ value (2*i+3)
            ulong[] bits = new ulong[(m + 63) >> 6];

            int sqrt = (int)Math.Sqrt((double)n);
            int maxI = (sqrt - 1) >> 1; // (2*maxI+3) ≤ sqrt(n)

            for (int i = 0; i <= maxI; i++)
            {
                if ((bits[i >> 6] & (1UL << (i & 63))) == 0)
                {
                    int p = (i << 1) + 3;
                    long start = 1L * p * p;         // first composite
                    int j = (int)((start - 3) >> 1); // index for odd number (2*j+3)
                    for (; j < m; j += p)
                        bits[j >> 6] |= 1UL << (j & 63);
                }
            }

            var res = new List<int>(m / 2 + 1) { 2 };

            for (int i = 0; i < m; i++)
                if ((bits[i >> 6] & (1UL << (i & 63))) == 0)
                    res.Add((i << 1) + 3);

            return res.ToArray();
        }

        /// <summary>
        /// Returns the radical of an integer.
        /// </summary>
        /// <param name="n">Value.</param>
        /// <returns>Integer number.</returns>
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
        /// <param name="n">Value.</param>
        /// <returns>Integer number.</returns>
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
        /// </summary>
        /// <remarks>
        /// Least-significant digit first: 10 in base 2 is {0,1,0,1}. Negative values are not supported.
        /// </remarks>
        /// <param name="x">Byte.</param>
        /// <param name="newbase">Base.</param>
        /// <returns>Array.</returns>
        public static int[] Decimal2Base(long x, int newbase)
        {
            if (x < 0) throw new ArgumentOutOfRangeException(nameof(x), "Digit arrays represent nonnegative integers.");
            if (newbase < 2) throw new ArgumentOutOfRangeException(nameof(newbase));
            var digits = new List<int>();
            do
            {
                digits.Add((int)(x % newbase));
                x /= newbase;
            } while (x != 0);
            return digits.ToArray();
        }
        /// <summary>
        /// Returns the decimal Number represented in decimal notation.
        /// </summary>
        /// <remarks>
        /// Least-significant digit first: {0,1,0,1} in base 2 represents 10. Overflow throws OverflowException.
        /// </remarks>
        /// <param name="x">Array.</param>
        /// <param name="thisbase">Base.</param>
        /// <returns>Integer number.</returns>
        public static long Base2Decimal(int[] x, int thisbase)
        {
            return AccumulateDigits(x, thisbase, true);
        }
        /// <summary>
        /// Returns a number that interprets the specified vector in decimal.
        /// </summary>
        /// <remarks>
        /// Most-significant digit first: {1,0,1,0} represents 1010. Overflow throws OverflowException.
        /// </remarks>
        /// <param name="x">Array.</param>
        /// <returns>Integer number.</returns>
        public static long Vector2Numeral(int[] x)
        {
            return AccumulateDigits(x, base10, false);
        }
        /// <summary>
        /// Returns a vector representing the decomposition of a decimal number into components.
        /// </summary>
        /// <remarks>
        /// Most-significant digit first: 1010 becomes {1,0,1,0}. Negative values are not supported.
        /// </remarks>
        /// <param name="x">Value.</param>
        /// <returns>Array.</returns>
        public static int[] Numeral2Vector(long x)
        {
            int[] digits = Decimal2Base(x, base10);
            Array.Reverse(digits);
            return digits;
        }
        /// <summary>
        /// Returns the digit count of the magnitude of an integer; zero has one digit.
        /// </summary>
        /// <param name="x">Byte.</param>
        /// <param name="numbase">Base.</param>
        /// <returns>Integer number.</returns>
        public static int NumLength(long x, int numbase)
        {
            if (numbase < 2) throw new ArgumentOutOfRangeException(nameof(numbase));
            ulong value = UnsignedMagnitude(x);
            int length = 1;
            while (value >= (ulong)numbase) { value /= (ulong)numbase; length++; }
            return length;
        }
        #endregion

        #region Solutions
        /// <summary>
        /// Returns the value of the hypotenuse.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Value.</returns>
        public static float Hypotenuse(float a, float b)
        {
            float r = 0.0f;
            float absA = Math.Abs(a);
            float absB = Math.Abs(b);

            if (absA > absB)
            {
                r = b / a;
                r = absA * Maths.Sqrt(1 + r * r);
            }
            else if (b != 0)
            {
                r = a / b;
                r = absB * Maths.Sqrt(1 + r * r);
            }

            return r;
        }
        /// <summary>
        /// Returns the value of the hypotenuse.
        /// </summary>
        /// <param name="z">Value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Hypotenuse(Complex32 z) => Maths.Hypotenuse(z.Real, z.Imag);

        /// <summary>
        /// Implements the solution of a cubic equation of the form:
        /// x^3 + a*x^2 + b*x + c = 0.
        /// </summary>
        /// <param name="a">Coefficient "a".</param>
        /// <param name="b">Coefficient "b".</param>
        /// <param name="c">Coefficient "c".</param>
        /// <returns>Array.</returns>
        public static Complex32[] Cubic(float a, float b, float c)
        {
            if (c == 0)
            {
                Complex32[] pair = Quadratic(1, a, b);
                return new[] { new Complex32(0, 0), pair[0], pair[1] };
            }
            if ((double)a * a > 1e8 * Math.Abs(b) && Math.Abs((double)a * a * a) > 1e8 * Math.Abs(c))
                return CubicWithSeparatedRoot(a, b, c);
            double shift = a / 3.0;
            double p = b - (double)a * a / 3.0;
            double q = 2.0 * a * a * a / 27.0 - (double)a * b / 3.0 + c;
            double halfQ = q / 2.0;
            double thirdP = p / 3.0;
            double discriminant = halfQ * halfQ + thirdP * thirdP * thirdP;
            if (discriminant < 0)
            {
                double radius = 2 * Math.Sqrt(-thirdP);
                double cosine = -halfQ / Math.Sqrt(-thirdP * thirdP * thirdP);
                double angle = Math.Acos(Math.Max(-1, Math.Min(1, cosine))) / 3.0;
                return new[] {
                    new Complex32((float)(radius * Math.Cos(angle) - shift), 0),
                    new Complex32((float)(radius * Math.Cos(angle + 2 * Math.PI / 3) - shift), 0),
                    new Complex32((float)(radius * Math.Cos(angle - 2 * Math.PI / 3) - shift), 0) };
            }
            // Choose the larger Cardano radicand and obtain the other term from u*v=-p/3.
            double u = RealCubeRoot(-halfQ - (halfQ < 0 ? -1 : 1) * Math.Sqrt(discriminant));
            double v = u == 0 ? 0 : -thirdP / u;
            double sum = u + v;
            double realRoot = sum - shift;
            // Recover a small real root lost by cancellation in Cardano's sum.
            for (int i = 0; i < 3; i++)
            {
                double derivative = (3 * realRoot + 2 * a) * realRoot + b;
                if (derivative == 0) break;
                realRoot -= (((realRoot + a) * realRoot + b) * realRoot + c) / derivative;
            }
            double imaginary = Math.Sqrt(3) * (u - v) / 2;
            return new[] {
                new Complex32((float)realRoot, 0),
                new Complex32((float)(-(a + realRoot) / 2), (float)imaginary),
                new Complex32((float)(-(a + realRoot) / 2), (float)-imaginary) };
        }
        /// <summary>
        /// Implements a solution to a quadratic equation of the form: 
        /// a*x^2 + b*x + c = 0.
        /// </summary>
        /// <param name="a">Coefficient "a".</param>
        /// <param name="b">Coefficient "b".</param>
        /// <param name="c">Coefficient "c".</param>
        /// <returns>Array.</returns>
        public static Complex32[] Quadratic(float a, float b, float c)
        {
            if (a == 0) throw new ArgumentOutOfRangeException(nameof(a), "The quadratic coefficient must be nonzero.");
            double discriminant = (double)b * b - 4.0 * a * c;
            if (discriminant < 0)
            {
                double re = -(double)b / (2.0 * a);
                double im = Math.Sqrt(-discriminant) / (2.0 * a);
                return new[] { (Complex32)new Complex(re, im), (Complex32)new Complex(re, -im) };
            }
            double q = -0.5 * (b + (b < 0 ? -1 : 1) * Math.Sqrt(discriminant));
            if (q == 0) return new[] { new Complex32(0, 0), new Complex32(0, 0) };
            return new[] { new Complex32((float)(q / a), 0), new Complex32((float)(c / q), 0) };
        }
        /// <summary>
        /// Implements the solution of a biquadratic equation of the form:
        /// a*x^4 + b*x^2 + c = 0.
        /// </summary>
        /// <param name="a">Coefficient "a".</param>
        /// <param name="b">Coefficient "b".</param>
        /// <param name="c">Coefficient "c".</param>
        /// <returns>Array.</returns>
        public static Complex32[] BiQuadratic(float a, float b, float c)
        {
            var s = Quadratic(a, b, c);
            return new Complex32[] {     Maths.Sqrt(s[0]),
                                      -Maths.Sqrt(s[0]),
                                       Maths.Sqrt(s[1]),
                                      -Maths.Sqrt(s[1]) };
        }
        #endregion

        #region Givens rotation
        /// <summary>
        /// Implements the construction of the Givens rotation matrix for a pair of real numbers.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Matrix.</returns>
        public static float[,] Rotation(float a, float b)
        {
            // MATLAB version of
            // Givens rotations:
            float c, s;
            float absx = Maths.Abs(a);

            if (absx == 0)
            {
                c = 0.0f;
                s = 1.0f;
            }
            else
            {
                float[] v = new float[] { a, b };
                float norm = v.Norm();
                c = absx / norm;
                s = a / absx * (b / norm);
            }

            return new float[,] { { c, s }, { -s, c } };
        }
        /// <summary>
        /// Implements the construction of the Givens rotation matrix for a pair of real numbers.
        /// </summary>
        /// <param name="a">Value.</param>
        /// <param name="b">Value.</param>
        /// <returns>Matrix.</returns>
        public static Complex32[,] Rotation(Complex32 a, Complex32 b)
        {
            // MATLAB version of
            // Givens rotations:
            Complex32 c, s;
            Complex32 absx = Maths.Abs(a);

            if (absx == 0)
            {
                c = 0.0;
                s = 1.0;
            }
            else
            {
                Complex32[] v = new Complex32[] { a, b };
                float norm = v.Norm();
                c = absx / norm;
                s = a / absx * (b.Conjugate / norm);
            }

            return new Complex32[,] { { c, s }, { -s.Conjugate, c } };
        }
        #endregion

        #region Other
        /// <summary>
        /// Returns the value of <c>|a|</c> with the sign of <paramref name="sign"/> (copysign).
        /// </summary>
        /// <param name="magnitude">Value providing the magnitude.</param>
        /// <param name="sign">Value providing the sign.</param>
        /// <returns>Value.</returns>
        public static float Sign(float magnitude, float sign)
        {
            return (sign >= 0.0) ? Math.Abs(magnitude) : -Math.Abs(magnitude);
        }
        /// <summary>
        /// Complex signum: returns z / |z| (unit complex) or 0 for z == 0.
        /// </summary>
        /// <param name="z">Complex value.</param>
        /// <returns>Value.</returns>
        public static Complex32 Sign(Complex32 z)
        {
            float re = z.Real, im = z.Imag;
            float r = Maths.Sqrt(re * re + im * im);
            if (r == 0f) return new Complex32(0f, 0f);
            return new Complex32(re / r, im / r);
        }

        /// <summary>
        /// Copies sign.
        /// </summary>
        /// <param name="magnitude">Value providing the magnitude.</param>
        /// <param name="sign">Value providing the sign.</param>
        /// <returns>Value.</returns>
        public static float CopySign(float magnitude, float sign)
        {
            return Math.Abs(magnitude) * Math.Sign(sign);
        }
        /// <summary>
        /// Copy phase from 'sign' to a real magnitude |magnitude|.
        /// If sign == 0, returns +|magnitude| on the real axis.
        /// </summary>
        /// <param name="magnitude">Value providing the magnitude.</param>
        /// <param name="sign">Value providing the sign.</param>
        /// <returns>Value.</returns>
        public static Complex32 CopySign(float magnitude, Complex32 sign)
        {
            float mag = Math.Abs(magnitude);
            float re = sign.Real, im = sign.Imag;
            float r = Maths.Sqrt(re * re + im * im);
            if (r == 0f) return new Complex32(mag, 0f);
            float s = mag / r;
            return new Complex32(re * s, im * s);
        }
        /// <summary>
        /// Copy phase from 'sign' to the magnitude |magnitude| of a complex number.
        /// If sign == 0, returns +|magnitude| on the real axis.
        /// </summary>
        /// <param name="magnitude">Value providing the magnitude.</param>
        /// <param name="sign">Value providing the sign.</param>
        /// <returns>Value.</returns>
        public static Complex32 CopySign(Complex32 magnitude, Complex32 sign)
        {
            float mag = Maths.Sqrt(magnitude.Real * magnitude.Real + magnitude.Imag * magnitude.Imag);
            float re = sign.Real, im = sign.Imag;
            float r = Maths.Sqrt(re * re + im * im);
            if (r == 0f) return new Complex32(mag, 0f);
            float s = mag / r;
            return new Complex32(re * s, im * s);
        }
        /// <summary>
        /// Copy phase from 'sign' to the magnitude |magnitude| of a complex number.
        /// If sign == 0, returns +|magnitude| on the real axis.
        /// </summary>
        /// <param name="magnitude">Value providing the magnitude.</param>
        /// <param name="sign">Value providing the sign.</param>
        /// <returns>Value.</returns>
        public static Complex32 CopySign(Complex32 magnitude, float sign)
        {
            if (float.IsNaN(sign))
                return Complex32.NaN;

            int s = Math.Sign(sign);
            if (s > 0) return magnitude;
            if (s < 0) return -magnitude;
            return Complex32.Zero;
        }

        /// <summary>
        /// Normalizes a variable relative to the {min, max} range.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="min">Minimum value.</param>
        /// <param name="max">Maximum value.</param>
        /// <returns>Value.</returns>
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
        /// <param name="x">Value.</param>
        /// <param name="min">Minimum value.</param>
        /// <param name="max">Maximum value.</param>
        /// <returns>Value.</returns>
        public static float Normalize(float x, float min, float max)
        {
            float a = max - min;
            float b = x - min;
            float c = (a != 0) ? b / a : x;
            return c;
        }
        #endregion

        #region Private arithmetic and number-theory helpers
        /// <summary>
        /// Evaluates complex sinh using double-precision real trigonometric factors.
        /// </summary>
        /// <param name="value">Input value.</param>
        /// <returns>The hyperbolic sine in double-precision complex arithmetic.</returns>
        private static Complex ComplexSinh(Complex32 value)
        {
            return new Complex(Math.Sinh(value.Real) * Math.Cos(value.Imag), Math.Cosh(value.Real) * Math.Sin(value.Imag));
        }

        /// <summary>
        /// Evaluates complex cosh using double-precision real trigonometric factors.
        /// </summary>
        /// <param name="value">Input value.</param>
        /// <returns>The hyperbolic cosine in double-precision complex arithmetic.</returns>
        private static Complex ComplexCosh(Complex32 value)
        {
            return new Complex(Math.Cosh(value.Real) * Math.Cos(value.Imag), Math.Sinh(value.Real) * Math.Sin(value.Imag));
        }

        /// <summary>
        /// Evaluates real asinh through a positive magnitude to avoid cancellation for negative arguments.
        /// </summary>
        /// <param name="value">Input value.</param>
        /// <returns>The real inverse hyperbolic sine, preserving the sign of the input.</returns>
        private static double RealAsinh(double value)
        {
            double magnitude = Math.Abs(value);
            if (magnitude < 1e-4 || double.IsInfinity(value)) return value;
            double result = Math.Log(magnitude + Math.Sqrt(magnitude * magnitude + 1));
            return value < 0 ? -result : result;
        }

        /// <summary>
        /// Evaluates real acosh for arguments at least one; smaller arguments return NaN.
        /// </summary>
        /// <param name="value">Real argument, at least one.</param>
        /// <returns>The nonnegative inverse hyperbolic cosine, or NaN below one.</returns>
        private static double RealAcosh(double value)
        {
            if (value < 1) return double.NaN;
            return Math.Log(value + Math.Sqrt((value - 1) * (value + 1)));
        }

        /// <summary>
        /// Evaluates log(1 + value) with a correction for rounding in the addition near zero.
        /// </summary>
        /// <param name="value">Finite real argument greater than -1; -1 gives the logarithmic pole.</param>
        /// <returns>The natural logarithm of 1 + value.</returns>
        private static double LogOnePlus(double value)
        {
            double sum = 1 + value;
            if (sum == 1) return value;
            return Math.Log(sum) * value / (sum - 1);
        }

        /// <summary>
        /// Evaluates the principal log(1 + value) using a corrected modulus and atan2 phase.
        /// </summary>
        /// <param name="value">Input value.</param>
        /// <returns>The principal complex logarithm of 1 + value.</returns>
        private static Complex ComplexLogOnePlus(Complex value)
        {
            double x = value.Real, y = value.Imaginary;
            return new Complex(0.5 * LogOnePlus(2 * x + x * x + y * y), Math.Atan2(y, 1 + x));
        }

        /// <summary>
        /// Evaluates complex atan with the branch convention used by the reciprocal arccotangent implementation.
        /// </summary>
        /// <remarks>On imaginary cuts with absolute imaginary part greater than one, the real part has the sign of the input imaginary component.</remarks>
        /// <param name="value">Input value.</param>
        /// <returns>The complex inverse tangent with the stated cut convention.</returns>
        private static Complex PrincipalAtan(Complex value)
        {
            double x = value.Real, y = value.Imaginary;
            double denominator = x * x + (y - 1) * (y - 1);
            double ratio = 4 * y / denominator;
            double imaginary = Math.Abs(ratio) < 0.5 ? LogOnePlus(ratio) :
                Math.Log(x * x + (y + 1) * (y + 1)) - Math.Log(denominator);
            double real = x == 0 && Math.Abs(y) > 1 ? Math.Sign(y) * Math.PI / 2 :
                0.5 * Math.Atan2(2 * x, 1 - x * x - y * y);
            return new Complex(real, 0.25 * imaginary);
        }

        /// <summary>
        /// Solves x^3 + a*x^2 + b*x + c = 0 when one root is isolated by a dominant quadratic coefficient.
        /// </summary>
        /// <remarks>Refines the isolated root and recovers the smaller roots from their sum and product, avoiding cancellation in the depressed cubic.</remarks>
        /// <param name="a">Coefficient of x^2 in the monic cubic.</param>
        /// <param name="b">Coefficient of x.</param>
        /// <param name="c">Constant coefficient.</param>
        /// <returns>The isolated real root and the two roots recovered from its residual quadratic.</returns>
        private static Complex32[] CubicWithSeparatedRoot(double a, double b, double c)
        {
            // When |a| dominates, depressing the cubic loses the two smaller roots.
            // Refine the isolated root and recover the remaining sum/product without cancellation.
            double root = -a;
            for (int i = 0; i < 4; i++)
                root -= (((root + a) * root + b) * root + c) / ((3 * root + 2 * a) * root + b);
            double product = -c / root;
            double sum = (b - product) / root;
            double discriminant = sum * sum - 4 * product;
            Complex first, second;
            if (discriminant < 0)
            {
                first = new Complex(sum / 2, Math.Sqrt(-discriminant) / 2);
                second = Complex.Conjugate(first);
            }
            else
            {
                double q = 0.5 * (sum + (sum < 0 ? -1 : 1) * Math.Sqrt(discriminant));
                first = q;
                second = q == 0 ? 0 : product / q;
            }
            return new[] { new Complex32((float)root, 0), (Complex32)first, (Complex32)second };
        }

        /// <summary>
        /// Evaluates real atanh using its small-argument limit and a logarithmic ratio.
        /// </summary>
        /// <param name="value">Real argument in [-1, 1].</param>
        /// <returns>The real inverse hyperbolic tangent, with infinite endpoint limits and NaN outside [-1, 1].</returns>
        private static double RealAtanh(double value)
        {
            if (Math.Abs(value) < 1e-4) return value;
            return 0.5 * Math.Log((1 + value) / (1 - value));
        }

        /// <summary>
        /// Returns the real cube root, preserving the sign of a negative argument.
        /// </summary>
        /// <param name="value">Input value.</param>
        /// <returns>The real cube root of value.</returns>
        private static double RealCubeRoot(double value)
        {
            double root = Math.Pow(Math.Abs(value), 1.0 / 3);
            return value < 0 ? -root : root;
        }

        /// <summary>
        /// Returns the unsigned magnitude of a signed integer, including Int64.MinValue.
        /// </summary>
        /// <param name="value">Input value.</param>
        /// <returns>The exact magnitude as an unsigned integer.</returns>
        private static ulong UnsignedMagnitude(long value)
        {
            return value < 0 ? (ulong)(-(value + 1)) + 1 : (ulong)value;
        }

        /// <summary>
        /// Computes a nonnegative greatest common divisor with the Euclidean remainder algorithm.
        /// </summary>
        /// <param name="a">First nonnegative integer.</param>
        /// <param name="b">Second nonnegative integer.</param>
        /// <returns>The nonnegative GCD; zero when both inputs are zero.</returns>
        private static ulong UnsignedGcd(ulong a, ulong b)
        {
            while (b != 0)
            {
                ulong remainder = a % b;
                a = b;
                b = remainder;
            }
            return a;
        }

        /// <summary>
        /// Computes a nonnegative GCD and Bezout coefficients using exact BigInteger arithmetic.
        /// </summary>
        /// <param name="a">First signed integer.</param>
        /// <param name="b">Second signed integer.</param>
        /// <returns>An array [gcd, x, y] satisfying a*x + b*y = gcd.</returns>
        private static BigInteger[] ExtendedGcd(BigInteger a, BigInteger b)
        {
            BigInteger oldR = a, r = b, oldS = 1, s = 0, oldT = 0, t = 1;
            while (r != 0)
            {
                BigInteger quotient = oldR / r;
                BigInteger nextR = oldR - quotient * r;
                BigInteger nextS = oldS - quotient * s;
                BigInteger nextT = oldT - quotient * t;
                oldR = r; r = nextR;
                oldS = s; s = nextS;
                oldT = t; t = nextT;
            }
            if (oldR < 0) { oldR = -oldR; oldS = -oldS; oldT = -oldT; }
            return new[] { oldR, oldS, oldT };
        }

        /// <summary>
        /// Multiplies two unsigned integers modulo a positive modulus without overflowing the product.
        /// </summary>
        /// <param name="a">First unsigned factor.</param>
        /// <param name="b">Second unsigned factor.</param>
        /// <param name="modulus">Positive modulus.</param>
        /// <returns>The remainder of the exact product a*b modulo modulus.</returns>
        private static ulong MultiplyModulo(ulong a, ulong b, ulong modulus)
        {
            if (a <= uint.MaxValue && b <= uint.MaxValue) return a * b % modulus;
            return (ulong)((BigInteger)a * b % modulus);
        }

        /// <summary>
        /// Evaluates an unsigned modular power by repeated squaring with exact intermediate products.
        /// </summary>
        /// <param name="value">Unsigned base.</param>
        /// <param name="exponent">Nonnegative integer exponent.</param>
        /// <param name="modulus">Modulus, greater than one at the call sites.</param>
        /// <returns>The modular power; exponent zero returns one.</returns>
        private static ulong PowerModulo(ulong value, ulong exponent, ulong modulus)
        {
            ulong result = 1;
            while (exponent != 0)
            {
                if ((exponent & 1) != 0) result = MultiplyModulo(result, value, modulus);
                exponent >>= 1;
                if (exponent != 0) value = MultiplyModulo(value, value, modulus);
            }
            return result;
        }

        // Sorenson and Webster, https://arxiv.org/abs/1509.00864:
        // the first composite passing all twelve bases exceeds the entire Int64 range.
        /// <summary>
        /// Miller–Rabin witnesses sufficient for primality decisions over the signed 64-bit domain.
        /// </summary>
        private static readonly uint[] PrimalityBases = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };

        /// <summary>
        /// Tests primality with Miller–Rabin bases sufficient throughout the signed 64-bit input domain.
        /// </summary>
        /// <param name="value">Nonnegative integer from the signed 64-bit input domain.</param>
        /// <returns>True for a prime; false for a composite integer, zero, or one.</returns>
        private static bool IsPrimeUnsigned(ulong value)
        {
            if (value < 2) return false;
            foreach (uint prime in PrimalityBases)
            {
                if (value % prime == 0) return value == prime;
            }
            ulong oddPart = value - 1;
            int powersOfTwo = 0;
            while ((oddPart & 1) == 0) { oddPart >>= 1; powersOfTwo++; }
            foreach (uint witness in PrimalityBases)
            {
                ulong residue = PowerModulo(witness, oddPart, value);
                if (residue == 1 || residue == value - 1) continue;
                bool passed = false;
                for (int i = 1; i < powersOfTwo; i++)
                {
                    residue = MultiplyModulo(residue, residue, value);
                    if (residue == value - 1) { passed = true; break; }
                }
                if (!passed) return false;
            }
            return true;
        }

        /// <summary>
        /// Appends the prime factors of a positive integer, retaining their multiplicities.
        /// </summary>
        /// <param name="value">Positive integer to factor.</param>
        /// <param name="factors">Destination list for prime factors, including repetitions.</param>
        private static void FactorInteger(ulong value, List<ulong> factors)
        {
            if (value == 1) return;
            if (IsPrimeUnsigned(value)) { factors.Add(value); return; }
            ulong divisor = FindDivisor(value);
            FactorInteger(divisor, factors);
            FactorInteger(value / divisor, factors);
        }

        /// <summary>
        /// Finds a proper divisor of a composite integer using Pollard–Brent retries and a trial-division fallback.
        /// </summary>
        /// <param name="value">Composite integer in the signed 64-bit input domain.</param>
        /// <returns>A divisor strictly between one and the composite input.</returns>
        private static ulong FindDivisor(ulong value)
        {
            foreach (uint prime in PrimalityBases)
                if (value % prime == 0) return prime;

            // Brent's batched Pollard rho. Restart failed cycles; a failed split is not primality evidence.
            // Bounded attempts plus exact trial division guarantee a finite fallback.
            for (ulong constant = 1; constant <= 32; constant++)
            {
                ulong y = 2, x = 0, saved = 0, divisor = 1;
                for (int length = 1; length <= 131072 && divisor == 1; length *= 2)
                {
                    x = y;
                    for (int i = 0; i < length; i++) y = (MultiplyModulo(y, y, value) + constant) % value;
                    for (int start = 0; start < length && divisor == 1; start += 64)
                    {
                        saved = y;
                        ulong product = 1;
                        int count = Math.Min(64, length - start);
                        for (int i = 0; i < count; i++)
                        {
                            y = (MultiplyModulo(y, y, value) + constant) % value;
                            ulong difference = x > y ? x - y : y - x;
                            product = MultiplyModulo(product, difference, value);
                        }
                        divisor = UnsignedGcd(product, value);
                    }
                }
                if (divisor == value)
                {
                    for (int i = 0; i < 64; i++)
                    {
                        saved = (MultiplyModulo(saved, saved, value) + constant) % value;
                        divisor = UnsignedGcd(x > saved ? x - saved : saved - x, value);
                        if (divisor != 1) break;
                    }
                }
                if (divisor > 1 && divisor < value) return divisor;
            }
            for (ulong divisor = 41; divisor <= value / divisor; divisor += 2)
                if (value % divisor == 0) return divisor;
            return value;
        }

        /// <summary>
        /// Decodes positional digits with checked accumulation and explicit digit-order selection.
        /// </summary>
        /// <param name="digits">Positional digits; each must be nonnegative and less than the radix.</param>
        /// <param name="radix">Integer radix, at least two.</param>
        /// <param name="leastSignificantFirst">True for least-significant-first digit order; false for most-significant-first order.</param>
        /// <returns>The decoded nonnegative Int64 value.</returns>
        private static long AccumulateDigits(int[] digits, int radix, bool leastSignificantFirst)
        {
            if (digits == null) throw new ArgumentNullException(nameof(digits));
            if (radix < 2) throw new ArgumentOutOfRangeException(nameof(radix));
            long result = 0;
            for (int i = 0; i < digits.Length; i++)
            {
                int digit = digits[leastSignificantFirst ? digits.Length - i - 1 : i];
                if (digit < 0 || digit >= radix) throw new ArgumentOutOfRangeException(nameof(digits), "Digits must be in the range [0, radix).");
                result = checked(result * radix + digit);
            }
            return result;
        }
        #endregion
    }
}
