using System;
using System.Collections.Generic;
using System.Linq;

namespace UMapx.Core
{
    /// <summary>
    /// Used to implement basic algebraic, trigonometric and hyperbolic operations.
    /// </summary>
    public static class Maths
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
        /// <param name="x">Value</param>
        /// <returns>Byte</returns>
        public static byte Byte(float x)
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
        /// Converts a value to an <see cref="sbyte"/> type and clamps it to the range [-128, 127].
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>sbyte</returns>
        public static sbyte sByte(float x)
        {
            return (sbyte)((x > 127) ? 127 : ((x < -128) ? -128 : x));
        }
        /// <summary>
        /// Converts a value to an <see cref="sbyte"/> type and clamps it to the range [-128, 127].
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>sbyte</returns>
        public static sbyte sByte(int x)
        {
            return (sbyte)((x > 127) ? 127 : ((x < -128) ? -128 : x));
        }

        /// <summary>
        /// Converts a value to a type float.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Float(float x)
        {
            return (x > 1.0f) ? 1.0f : ((x < 0) ? 0 : x);
        }
        /// <summary>
        /// Checks if value is in the specified range.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="xmax">Maximum value</param>
        /// <returns>Boolean</returns>
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
        /// <param name="a">Value</param>
        /// <returns>Boolean</returns>
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
        /// <param name="a">Complex number</param>
        /// <returns>Boolean</returns>
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
        /// <param name="n">Integer number</param>
        /// <returns>Boolean</returns>
        public static bool IsSquare(float n)
        {
            float sq = (int)Math.Sqrt(n);
            return (sq * sq == n);
        }
        /// <summary>
        /// Checks whether a number is a power of another number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Boolean</returns>
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
        /// <param name="a">Value</param>
        /// <returns>Boolean</returns>
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
        /// <param name="a">Value</param>
        /// <returns>Boolean</returns>
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
        /// <param name="a">Value</param>
        /// <returns>Boolean</returns>
        public static bool IsNotEven(float a)
        {
            return !IsEven(a);
        }
        /// <summary>
        /// Returns the number raised to the second power.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Pow(float a)
        {
            return (float)Math.Pow(a, 2);
        }
        /// <summary>
        /// Returns the number raised to the power.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Power</param>
        /// <returns>Value</returns>
        public static float Pow(float a, float b)
        {
            return (float)Math.Pow(a, b);
        }
        /// <summary>
        /// Returns the exponent raised to the power.
        /// </summary>
        /// <param name="a">Power</param>
        /// <returns>Value</returns>
        public static float Exp(float a)
        {
            return (float)Math.Pow(E, a);
        }
        /// <summary>
        /// Returns the natural logarithm of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Log(float a)
        {
            return (float)Math.Log(a);
        }
        /// <summary>
        /// Returns the decimal logarithm of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Log10(float a)
        {
            return (float)Math.Log(a, 10.0);
        }
        /// <summary>
        /// Returns the binary logarithm of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Log2(float a)
        {
            return (float)Math.Log(a, 2);
        }
        /// <summary>
        /// Returns the logarithm of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Base</param>
        /// <returns>Value</returns>
        public static float Log(float a, float b)
        {
            return (float)Math.Log(a, b);
        }
        /// <summary>
        /// Returns the square root of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Sqrt(float a)
        {
            return (float)Math.Sqrt(a);
        }
        /// <summary>
        /// Returns the root of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Power</param>
        /// <returns>Value</returns>
        public static float Sqrt(float a, float b)
        {
            return (float)Math.Pow(a, 1.0 / b);
        }
        /// <summary>
        /// Returns the modulus of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
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
        /// <param name="a">First number</param>
        /// <param name="b">Second number</param>
        /// <returns>Value</returns>
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
        /// <param name="a">First number</param>
        /// <param name="b">Second number</param>
        /// <param name="c">Third number</param>
        /// <returns>Value</returns>
        public static float Max(float a, float b, float c)
        {
            return Max(a, Max(b, c));
        }
        /// <summary>
        /// Returns the smallest of two numbers.
        /// </summary>
        /// <param name="a">First number</param>
        /// <param name="b">Second number</param>
        /// <returns>Value</returns>
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
        /// <param name="a">First number</param>
        /// <param name="b">Second number</param>
        /// <param name="c">Third number</param>
        /// <returns>Value</returns>
        public static float Min(float a, float b, float c)
        {
            return Min(a, Min(b, c));
        }
        /// <summary>
        /// Returns the sign of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
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
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Floor(float a)
        {
            return (float)Math.Floor(a);
        }
        /// <summary>
        /// Returns the rounded number up.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Ceil(float a)
        {
            return (float)Math.Ceiling(a);
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Round(float a)
        {
            return (float)Math.Round(a, 0);
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="dig">Digits</param>
        /// <returns>Value</returns>
        public static float Round(float a, int dig)
        {
            return (float)Math.Round(a, dig);
        }
        #endregion

        #region Complex number
        /// <summary>
        /// Returns the modulus of a complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static float Abs(Complex32 a)
        {
            return a.Abs;
        }
        /// <summary>
        /// Returns the angle of a complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static float Angle(Complex32 a)
        {
            return a.Angle;
        }
        /// <summary>
        /// Returns the natural logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Log(Complex32 a)
        {
            return Log(a, E);
        }
        /// <summary>
        /// Returns the decimal logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Log10(Complex32 a)
        {
            return Log(a, 10.0f);
        }
        /// <summary>
        /// Returns the binary logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Log2(Complex32 a)
        {
            return Log(a, 2.0f);
        }
        /// <summary>
        /// Returns the logarithm of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Base</param>
        /// <returns>Complex number</returns>
        public static Complex32 Log(Complex32 a, float b)
        {
            return new Complex32((float)Math.Log(a.Abs), a.Angle) / Math.Log(b);
        }
        /// <summary>
        /// Returns the exponent raised to a complex degree.
        /// </summary>
        /// <param name="a">Power</param>
        /// <returns>Complex number</returns>
        public static Complex32 Exp(Complex32 a)
        {
            return Pow(E, a);
        }
        /// <summary>
        /// Returns the number raised to a complex power.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Power</param>
        /// <returns>Complex number</returns>
        public static Complex32 Pow(float a, Complex32 b)
        {
            float r = (float)Math.Pow(a, b.Real);
            return new Complex32(r * (float)Math.Cos(b.Imag), r * (float)Math.Sin(b.Imag));
        }
        /// <summary>
        /// Returns the number raised to the power.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Power</param>
        /// <returns>Complex number</returns>
        public static Complex32 Pow(Complex32 a, float b)
        {
            return Math.Pow(a.Abs, b) * new Complex32((float)Math.Cos(b * a.Angle), (float)Math.Sin(b * a.Angle));
        }
        /// <summary>
        /// Returns the number raised to the power.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Power</param>
        /// <returns>Complex number</returns>
        public static Complex32 Pow(Complex32 a, Complex32 b)
        {
            return Maths.Exp(b * Maths.Log(a));
        }
        /// <summary>
        /// Returns the square root of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Sqrt(Complex32 a)
        {
            return Maths.Sqrt(a, 2);
        }
        /// <summary>
        /// Returns the root of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Power</param>
        /// <returns>Complex number</returns>
        public static Complex32 Sqrt(Complex32 a, float b)
        {
            return Maths.FromPolar((float)Math.Sqrt(a.Abs), a.Angle / b);
        }
        /// <summary>
        /// Returns the root of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Power</param>
        /// <returns>Complex number</returns>
        public static Complex32 Sqrt(Complex32 a, Complex32 b)
        {
            return Maths.Exp(Maths.Log(a) / b);
        }
        /// <summary>
        /// Returns complex number.
        /// </summary>
        /// <param name="abs">Modulus</param>
        /// <param name="angle">Angle</param>
        /// <returns>Complex number</returns>
        public static Complex32 FromPolar(float abs, float angle)
        {
            return new Complex32(abs * (float)Math.Cos(angle), abs * (float)Math.Sin(angle));
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Round(Complex32 a)
        {
            return Maths.Round(a, 0);
        }
        /// <summary>
        /// Returns the rounded number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="dig">Digits</param>
        /// <returns>Complex number</returns>
        public static Complex32 Round(Complex32 a, int dig)
        {
            return new Complex32((float)Math.Round(a.Real, dig), (float)Math.Round(a.Imag, dig));
        }
        #endregion
        #endregion

        #region Trigonometric
        #region Real number
        /// <summary>
        /// Returns the cosine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Cos(float a)
        {
            return (float)Math.Cos(a);
        }
        /// <summary>
        /// Returns the sine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Sin(float a)
        {
            return (float)Math.Sin(a);
        }
        /// <summary>
        /// Returns the tangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Tan(float a)
        {
            return (float)Math.Sin(a) / (float)Math.Cos(a);
        }
        /// <summary>
        /// Returns the cotangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Ctan(float a)
        {
            return (float)Math.Cos(a) / (float)Math.Sin(a);
        }
        /// <summary>
        /// Returns the secant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Sec(float a)
        {
            return 1.0f / (float)Math.Cos(a);
        }
        /// <summary>
        /// Returns the cosecant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Cosc(float a)
        {
            return 1.0f / (float)Math.Sin(a);
        }
        /// <summary>
        /// Returns the arcsine of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Asin(float a)
        {
            return (float)Math.Asin(a);
        }
        /// <summary>
        /// Returns the arccosine of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Acos(float a)
        {
            return (float)Math.Acos(a);
        }
        /// <summary>
        /// Returns the arctangent of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Atan(float a)
        {
            return (float)Math.Atan(a);
        }
        /// <summary>
        /// Returns the arccotangent of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Actan(float a)
        {
            return Pi / 2 - (float)Math.Atan(a);
        }
        /// <summary>
        /// Returns the arcsecance of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Asec(float a)
        {
            return (float)Math.Acos(1.0 / a);
        }
        /// <summary>
        /// Returns the arccosecant of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Acosc(float a)
        {
            return (float)Math.Asin(1.0 / a);
        }
        #endregion

        #region Complex number
        /// <summary>
        /// Returns the cosine of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Cos(Complex32 a)
        {
            return new Complex32((float)Math.Cos(a.Real) * (float)Math.Cosh(a.Imag), -((float)Math.Sin(a.Real) * (float)Math.Sinh(a.Imag)));
        }
        /// <summary>
        /// Returns the sine of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Sin(Complex32 a)
        {
            return new Complex32((float)Math.Sin(a.Real) * (float)Math.Cosh(a.Imag), (float)Math.Cos(a.Real) * (float)Math.Sinh(a.Imag));
        }
        /// <summary>
        /// Returns the tangent of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Tan(Complex32 a)
        {
            return Maths.Sin(a) / Maths.Cos(a);
        }
        /// <summary>
        /// Returns the cotangent of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Ctan(Complex32 a)
        {
            return Maths.Cos(a) / Maths.Sin(a);
        }
        /// <summary>
        /// Returns the secant of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Sec(Complex32 a)
        {
            return 1.0 / Maths.Cos(a);
        }
        /// <summary>
        /// Returns the cosecant of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Cosc(Complex32 a)
        {
            return 1.0 / Maths.Sin(a);
        }
        /// <summary>
        /// Returns the arccosine of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Acos(Complex32 a)
        {
            return -I * Maths.Log(a + I * Maths.Sqrt(1.0 - a * a));
        }
        /// <summary>
        /// Returns the arcsine of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Asin(Complex32 a)
        {
            return -I * Maths.Log(I * a + Maths.Sqrt(1.0 - a * a));
        }
        /// <summary>
        /// Returns the arctangent of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Atan(Complex32 a)
        {
            return I / 2.0 * (Maths.Log(1.0 - I * a) - Maths.Log(1.0 + I * a));
        }
        /// <summary>
        /// Returns the arccotangent of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Actan(Complex32 a)
        {
            return I / 2.0 * (Maths.Log((a - I) / a) - Maths.Log((a + I) / a));
        }
        /// <summary>
        /// Returns the arcsecance of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Asec(Complex32 a)
        {
            return Maths.Acos(1.0 / a);
        }
        /// <summary>
        /// Returns the arccosecant of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
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
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Sinh(float a)
        {
            return (float)Math.Sinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cosine of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Cosh(float a)
        {
            return (float)Math.Cosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic tangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Tanh(float a)
        {
            return (float)Math.Sinh(a) / (float)Math.Cosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cotangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Ctanh(float a)
        {
            return (float)Math.Cosh(a) / (float)Math.Sinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic secant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Sech(float a)
        {
            return 1.0f / (float)Math.Cosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cosecant of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Cosch(float a)
        {
            return 1.0f / (float)Math.Sinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic arcsine of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Asinh(float a)
        {
            return (float)Math.Log(a + Math.Sqrt(a * a + 1));
        }
        /// <summary>
        /// Returns the hyperbolic arccosine of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Acosh(float a)
        {
            if (a >= 0)
            {
                return (float)Math.Log(a + Math.Sqrt(a * a - 1));
            }
            return 0;
        }
        /// <summary>
        /// Returns the hyperbolic arctangent of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Atanh(float a)
        {
            return 1.0f / 2.0f * (float)Math.Log((1 + a) / (1 - a));
        }
        /// <summary>
        /// Returns the hyperbolic arccotangent of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Actanh(float a)
        {
            return 1.0f / 2.0f * (float)Math.Log((a + 1) / (a - 1));
        }
        /// <summary>
        /// Returns the hyperbolic arcsecance of a number.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Value</returns>
        public static float Asech(float a)
        {
            return (float)Math.Log((1 + (float)Math.Sqrt(1 - a * a)) / a);
        }
        /// <summary>
        /// Returns the hyperbolic arccosecant of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Acosch(float a)
        {
            if (a < 0)
            {
                return (float)Math.Log((1 - (float)Math.Sqrt(1 + a * a)) / a);
            }
            if (a > 0)
            {
                return (float)Math.Log((1 + (float)Math.Sqrt(1 + a * a)) / a);
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
        public static Complex32 Sinh(Complex32 a)
        {
            return new Complex32((float)Math.Sinh(a.Real) * (float)Math.Cos(a.Imag), (float)Math.Cosh(a.Real) * (float)Math.Sin(a.Imag));
        }
        /// <summary>
        /// Returns the hyperbolic cosine of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Cosh(Complex32 a)
        {
            return new Complex32((float)Math.Cosh(a.Real) * (float)Math.Cos(a.Imag), (float)Math.Sinh(a.Real) * (float)Math.Sin(a.Imag));
        }
        /// <summary>
        /// Returns the hyperbolic tangent of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Tanh(Complex32 a)
        {
            return Maths.Sinh(a) / Maths.Cosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cotangent of an angle.
        /// </summary>
        /// <param name="a">Angle in radians</param>
        /// <returns>Complex number</returns>
        public static Complex32 Ctanh(Complex32 a)
        {
            return Maths.Cosh(a) / Maths.Sinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic secant of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Sech(Complex32 a)
        {
            return 1.0 / Maths.Cosh(a);
        }
        /// <summary>
        /// Returns the hyperbolic cosecant of an angle.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Cosch(Complex32 a)
        {
            return 1.0 / Maths.Sinh(a);
        }
        /// <summary>
        /// Returns the hyperbolic arcsine of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Asinh(Complex32 a)
        {
            return Maths.Log(a + Maths.Sqrt(a * a + 1.0));
        }
        /// <summary>
        /// Returns the hyperbolic arccosine of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Acosh(Complex32 a)
        {
            return Maths.Log(a + Maths.Sqrt(a * a - 1.0));
        }
        /// <summary>
        /// Returns the hyperbolic arctangent of a number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <returns>Complex number</returns>
        public static Complex32 Atanh(Complex32 a)
        {
            return 1.0 / 2.0 * Maths.Log((1.0 + a) / (1.0 - a));
        }
        /// <summary>
        /// Returns the hyperbolic arccotangent of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Actanh(Complex32 a)
        {
            return 1.0 / 2.0 * Maths.Log((a + 1.0) / (a - 1.0));
        }
        /// <summary>
        /// Returns the hyperbolic arcsecance of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Asech(Complex32 a)
        {
            return Maths.Log(1.0 / a + Maths.Sqrt(1.0 / a + 1.0) + Maths.Sqrt(1.0 / a - 1.0));
        }
        /// <summary>
        /// Returns the hyperbolic arccosecant of a number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 Acosch(Complex32 a)
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
        /// <param name="p">Value</param>
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
        /// <param name="p">Value</param>
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
        /// <param name="a">Value</param>
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
        /// <param name="a">Value</param>
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
        /// <param name="a">Value</param>
        /// <param name="n">Modulo</param>
        /// <returns>Integer number</returns>
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
        /// <param name="a">Value</param>
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
        /// <param name="a">Value</param>
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
        /// Computes modular exponentiation using the left-to-right binary method.
        /// </summary>
        /// <param name="a">Base value</param>
        /// <param name="x">Exponent</param>
        /// <param name="p">Modulus</param>
        /// <returns>Result of a^x mod p</returns>
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
        /// Computes modular exponentiation using the right-to-left binary method.
        /// </summary>
        /// <param name="a">Base value</param>
        /// <param name="x">Exponent</param>
        /// <param name="p">Modulus</param>
        /// <returns>Result of a^x mod p</returns>
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
        /// <param name="a">Value</param>
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
        /// <param name="a">Value</param>
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
        /// <param name="a">Value</param>
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
        /// <param name="a">Value</param>
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
        /// <param name="n">Value</param>
        /// <returns>Value</returns>
        public static int Etf(int n)
        {
            // factorization with only primes
            int[] itf = Maths.Itf(n, true);
            float radical = 1;
            int length = itf.Length;

            // calculation radical
            for (int i = 0; i < length; i++)
            {
                radical *= 1.0f - 1.0f / itf[i];
            }
            return (int)(n * radical);
        }
        /// <summary>
        /// Returns the value of the Euler function.
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Value</returns>
        public static long Etf(long n)
        {
            // factorization with only primes
            long[] itf = Maths.Itf(n, true);
            float radical = 1;
            int length = itf.Length;

            // calculation radical
            for (int i = 0; i < length; i++)
            {
                radical *= 1.0f - 1.0f / itf[i];
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
        /// <param name="limit">Value</param>
        /// <returns>Array</returns>
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
            for (int low = 3; low <= limit;)
            {
                long highL = low + ((long)segmentOddCount << 1); // convert odd-count→integers
                if (highL > limitPlus1) highL = limitPlus1;
                int high = (int)highL; // exclusive

                // First odd ≥ low
                int firstOdd = (low | 1);

                // Number of odd integers in [firstOdd, high):
                // count = ceil((high - firstOdd)/2) = (high - firstOdd + 1) >> 1, clamped to ≥0
                int oddCount = high > firstOdd ? ((high - firstOdd + 1) >> 1) : 0;

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
                        int n = firstOdd + (i << 1);
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
        /// <param name="n">Value</param>
        /// <returns>Array</returns>
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
        /// <param name="n">Value</param>
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
        /// <param name="n">Value</param>
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
        /// <param name="x">Value</param>
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
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Value</returns>
        public static float Hypotenuse(float a, float b)
        {
            float r = 0.0f;
            float absA = Math.Abs(a);
            float absB = Math.Abs(b);

            if (absA > absB)
            {
                r = b / a;
                r = absA * (float)Math.Sqrt(1 + r * r);
            }
            else if (b != 0)
            {
                r = a / b;
                r = absB * (float)Math.Sqrt(1 + r * r);
            }

            return r;
        }
        /// <summary>
        /// Returns the value of the hypotenuse.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Value</returns>
        public static Complex32 Hypotenuse(Complex32 z) => Maths.Hypotenuse(z.Real, z.Imag);

        /// <summary>
        /// Implements the solution of a cubic equation of the form:
        /// x^3 + a*x^2 + b*x + c = 0.
        /// </summary>
        /// <param name="a">Coefficient "a"</param>
        /// <param name="b">Coefficient "b"</param>
        /// <param name="c">Coefficient "c"</param>
        /// <returns>Array</returns>
        public static Complex32[] Cubic(float a, float b, float c)
        {
            Complex32 x1 = 0, x2 = 0, x3 = 0;
            float Q = (a * a - 3.0f * b) / 9.0f;
            float R = (2.0f * a * a * a - 9.0f * a * b + 27.0f * c) / 54.0f;
            float S = Q * Q * Q - R * R;
            float a3 = a / 3.0f;
            float fi, v0, v1;

            if (S > 0)
            {
                fi = (float)Math.Acos(R / (float)Math.Pow(Q, 3.0f / 2.0f)) / 3.0f;
                v0 = -2 * (float)Math.Sqrt(Q);
                v1 = 2.0f / 3 * (float)Math.PI;

                x1 = v0 * Math.Cos(fi) - a3;
                x2 = v0 * Math.Cos(fi + v1) - a3;
                x3 = v0 * Math.Cos(fi - v1) - a3;
            }
            else if (S < 0)
            {
                if (Q > 0)
                {
                    fi = Maths.Acosh(Math.Abs(R) / (float)Math.Pow(Math.Abs(Q), 3.0f / 2.0f)) / 3.0f;
                    v0 = Math.Sign(R) * Maths.Sqrt(Q) * Maths.Cosh(fi);
                    v1 = Maths.Sqrt(3) * Maths.Sqrt(Q) * Maths.Sinh(fi);

                    x1 = -2 * v0 - a3;
                    x2 = v0 - a3 + Maths.I * v1;
                    x3 = x2.Conjugate;
                }
                else if (Q < 0)
                {
                    fi = Maths.Asinh(Math.Abs(R) / Maths.Pow(Math.Abs(Q), 3.0f / 2.0f)) / 3.0f;
                    v0 = Math.Sign(R) * Maths.Sqrt(Math.Abs(Q)) * Maths.Sinh(fi);
                    v1 = Maths.Sqrt(3) * Maths.Sqrt(Math.Abs(Q)) * Maths.Cosh(fi);

                    x1 = -2 * v0 - a3;
                    x2 = v0 - a3 + Maths.I * v1;
                    x3 = x2.Conjugate;
                }
                else if (Q == 0)
                {
                    x1 = -Maths.Sqrt(c - a * a * a / 27.0f, 3.0f) - a3;
                    v0 = Maths.Abs((a - 3 * x1) * (a + x1) - 4 * b);
                    x2 = Maths.I / 2.0 * Math.Sqrt(v0) - (a + x1) / 2.0;
                    x3 = x2.Conjugate;
                }
            }
            else if (S == 0)
            {
                v0 = (float)Math.Pow(R, 1.0 / 3.0);
                x1 = -2 * v0 - a3;
                x2 = x3 = v0 - a3;
            }
            return new Complex32[] { x1, x2, x3 };
        }
        /// <summary>
        /// Implements a solution to a quadratic equation of the form: 
        /// a*x^2 + b*x + c = 0.
        /// </summary>
        /// <param name="a">Coefficient "a"</param>
        /// <param name="b">Coefficient "b"</param>
        /// <param name="c">Coefficient "c"</param>
        /// <returns>Array</returns>
        public static Complex32[] Quadratic(float a, float b, float c)
        {
            float dis = b * b - 4 * a * c;
            float abs = (float)Math.Sqrt(Math.Abs(dis));
            Complex32 root = dis < 0 ? new Complex32(0, abs) : new Complex32(abs, 0);
            Complex32 q = -0.5 * (b + Math.Sign(b) * root);
            return new Complex32[] { q / a, c / q };
        }
        /// <summary>
        /// Implements the solution of a biquadratic equation of the form:
        /// a*x^4 + b*x^2 + c = 0.
        /// </summary>
        /// <param name="a">Coefficient "a"</param>
        /// <param name="b">Coefficient "b"</param>
        /// <param name="c">Coefficient "c"</param>
        /// <returns>Array</returns>
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
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Matrix</returns>
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
        /// <param name="a">Value</param>
        /// <param name="b">Value</param>
        /// <returns>Matrix</returns>
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
        /// <param name="magnitude">Value providing the magnitude</param>
        /// <param name="sign">Value providing the sign</param>
        /// <returns>Value</returns>
        public static float Sign(float magnitude, float sign)
        {
            return (sign >= 0.0) ? Math.Abs(magnitude) : -Math.Abs(magnitude);
        }
        /// <summary>
        /// Complex signum: returns z / |z| (unit complex) or 0 for z == 0.
        /// </summary>
        /// <param name="z">Complex value</param>
        /// <returns>Value</returns>
        public static Complex32 Sign(Complex32 z)
        {
            float re = z.Real, im = z.Imag;
            float r = (float)Math.Sqrt(re * re + im * im);
            if (r == 0f) return new Complex32(0f, 0f);
            return new Complex32(re / r, im / r);
        }

        /// <summary>
        /// Copies sign.
        /// </summary>
        /// <param name="magnitude">Value providing the magnitude</param>
        /// <param name="sign">Value providing the sign</param>
        /// <returns>Value</returns>
        public static float CopySign(float magnitude, float sign)
        {
            return Math.Abs(magnitude) * Math.Sign(sign);
        }
        /// <summary>
        /// Copy phase from 'sign' to a real magnitude |magnitude|.
        /// If sign == 0, returns +|magnitude| on the real axis.
        /// </summary>
        /// <param name="magnitude">Value providing the magnitude</param>
        /// <param name="sign">Value providing the sign</param>
        /// <returns>Value</returns>
        public static Complex32 CopySign(float magnitude, Complex32 sign)
        {
            float mag = Math.Abs(magnitude);
            float re = sign.Real, im = sign.Imag;
            float r = (float)Math.Sqrt(re * re + im * im);
            if (r == 0f) return new Complex32(mag, 0f);
            float s = mag / r;
            return new Complex32(re * s, im * s);
        }
        /// <summary>
        /// Copy phase from 'sign' to the magnitude |magnitude| of a complex number.
        /// If sign == 0, returns +|magnitude| on the real axis.
        /// </summary>
        /// <param name="magnitude">Value providing the magnitude</param>
        /// <param name="sign">Value providing the sign</param>
        /// <returns>Value</returns>
        public static Complex32 CopySign(Complex32 magnitude, Complex32 sign)
        {
            float mag = (float)Math.Sqrt(magnitude.Real * magnitude.Real + magnitude.Imag * magnitude.Imag);
            float re = sign.Real, im = sign.Imag;
            float r = (float)Math.Sqrt(re * re + im * im);
            if (r == 0f) return new Complex32(mag, 0f);
            float s = mag / r;
            return new Complex32(re * s, im * s);
        }
        /// <summary>
        /// Copy phase from 'sign' to the magnitude |magnitude| of a complex number.
        /// If sign == 0, returns +|magnitude| on the real axis.
        /// </summary>
        /// <param name="magnitude">Value providing the magnitude</param>
        /// <param name="sign">Value providing the sign</param>
        /// <returns>Value</returns>
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
        /// <param name="x">Value</param>
        /// <param name="min">Minimum value</param>
        /// <param name="max">Maximum value</param>
        /// <returns>Value</returns>
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
        /// <returns>Value</returns>
        public static float Normalize(float x, float min, float max)
        {
            float a = max - min;
            float b = x - min;
            float c = (a != 0) ? b / a : x;
            return c;
        }
        #endregion
    }
}
