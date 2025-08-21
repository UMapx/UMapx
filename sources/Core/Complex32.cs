using System;
using System.Numerics;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a complex number.
    /// </summary>
    [Serializable]
    public struct Complex32 : ICloneable
    {
        #region Private data
        /// <summary>
        /// The real part of the complex number.
        /// </summary>
        public float Real;
        /// <summary>
        /// The imaginary part of a complex number.
        /// </summary>
        public float Imag;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes the complex number.
        /// </summary>
        /// <param name="re">Real part of the complex number</param>
        /// <param name="im">Imaginary part of a complex number</param>
        public Complex32(float re, float im)
        {
            this.Real = re;
            this.Imag = im;
        }
        /// <summary>
        /// Gets the value of the module.
        /// </summary>
        public float Abs
        {
            get
            {
                return (float)Math.Sqrt(Real * Real + Imag * Imag);
            }
        }
        /// <summary>
        /// Gets the value of the module * module.
        /// </summary>
        public float Abs2
        {
            get
            {
                return Real * Real + Imag * Imag;
            }
        }
        /// <summary>
        /// Gets the value of the phase.
        /// </summary>
        public float Angle
        {
            get
            {
                return (float)Math.Atan2(Imag, Real);
            }
        }
        #endregion

        #region Overrides
        /// <summary>
        /// Returns the hash code for this object.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return this.Real.GetHashCode() ^ this.Imag.GetHashCode();
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the given value of type Complex.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is Complex32) ? (this == (Complex32)obj) : false;
        }
        /// <summary>
        /// Converts complex number to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return this.ToString("G6");
        }
        /// <summary>
        /// Converts complex number to its corresponding string representation.
        /// </summary>
        /// <param name="format">Format string</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public string ToString(string format)
        {
            return StringOptions.Disp(new float[] { this.Real, this.Imag }, format, StringOptions.C);
        }
        #endregion

        #region Bools
        /// <summary>
        /// Checks if two complex numbers are equal.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(Complex32 a, Complex32 b)
        {
            return ((a.Real == b.Real) && (a.Imag == b.Imag));
        }
        /// <summary>
        /// Checks if two complex numbers are not equal.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(Complex32 a, Complex32 b)
        {
            return !(a == b);
        }
        #endregion

        #region Operators
        /// <summary>
        /// The sum of two complex numbers.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator +(Complex32 a, Complex32 b)
        {
            return new Complex32(a.Real + b.Real, a.Imag + b.Imag);
        }
        /// <summary>
        /// The sum of a complex number and a real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Value</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator +(Complex32 a, float b)
        {
            return new Complex32(a.Real + b, a.Imag);
        }
        /// <summary>
        /// The sum of a complex number and a real number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator +(float a, Complex32 b)
        {
            return new Complex32(b.Real + a, b.Imag);
        }


        /// <summary>
        /// The difference of two complex numbers.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator -(Complex32 a, Complex32 b)
        {
            return new Complex32(a.Real - b.Real, a.Imag - b.Imag);
        }
        /// <summary>
        /// The difference between a complex number and a real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Value</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator -(Complex32 a, float b)
        {
            return new Complex32(a.Real - b, a.Imag);
        }
        /// <summary>
        /// The difference between a complex number and a real number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator -(float a, Complex32 b)
        {
            return new Complex32(a - b.Real, b.Imag);
        }
        /// <summary>
        /// Inverts complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator -(Complex32 a)
        {
            return new Complex32(-a.Real, -a.Imag);
        }


        /// <summary>
        /// Multiplies one complex number by another.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator *(Complex32 a, Complex32 b)
        {
            float aRe = a.Real, aIm = a.Imag;
            float bRe = b.Real, bIm = b.Imag;

            return new Complex32(aRe * bRe - aIm * bIm, aRe * bIm + aIm * bRe);
        }
        /// <summary>
        /// Multiplies real number by complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Value</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator *(float a, Complex32 b)
        {
            return new Complex32(b.Real * a, b.Imag * a);
        }
        /// <summary>
        /// Multiplies complex number by real number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator *(Complex32 a, float b)
        {
            return new Complex32(a.Real * b, a.Imag * b);
        }


        /// <summary>
        /// Divides one complex number by another.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator /(Complex32 a, Complex32 b)
        {
            float aRe = a.Real, aIm = a.Imag;
            float bRe = b.Real, bIm = b.Imag;
            float abs = bRe * bRe + bIm * bIm;
            float inv = 1 / abs;

            return new Complex32((aRe * bRe + aIm * bIm) * inv, (aIm * bRe - aRe * bIm) * inv);
        }
        /// <summary>
        /// Divides complex number by real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Value</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator /(Complex32 a, float b)
        {
            return new Complex32(a.Real / b, a.Imag / b);
        }
        /// <summary>
        /// Divides real number by complex number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex32 operator /(float a, Complex32 b)
        {
            // (a + 0i) / (bRe + i*bIm) = a*(bRe - i*bIm) / (bRe^2 + bIm^2)
            float bRe = b.Real;
            float bIm = b.Imag;
            float abs = bRe * bRe + bIm * bIm;

            return new Complex32(a * bRe / abs, -a * bIm / abs);
        }
        #endregion

        #region Conversion operators
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(Complex value)
        {
            return new Complex32((float)value.Real, (float)value.Imaginary);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(double value)
        {
            return new Complex32((float)value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(float value)
        {
            return new Complex32(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(long value)
        {
            return new Complex32(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(ulong value)
        {
            return new Complex32(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(short value)
        {
            return new Complex32(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(ushort value)
        {
            return new Complex32(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(int value)
        {
            return new Complex32(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(uint value)
        {
            return new Complex32(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(byte value)
        {
            return new Complex32(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(sbyte value)
        {
            return new Complex32(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex32(decimal value)
        {
            return new Complex32((float)value, 0);
        }
        #endregion

        #region Parsing
        /// <summary>
        /// Parses the string to complex number.
        /// <remarks>
        /// Example: "1 + 2i", "0.321 + 11i", ".1i".
        /// </remarks>
        /// </summary>
        /// <param name="s">Input string</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static Complex32 Parse(string s)
        {
            return StringOptions.Compar(s);
        }
        /// <summary>
        /// Tries to parse the string to complex number.
        /// </summary>
        /// <param name="complex">Input string</param>
        /// <param name="result">Complex number</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string complex, out Complex32 result)
        {
            try
            {
                result = Complex32.Parse(complex);
                return true;
            }
            catch (FormatException)
            {
                result = new Complex32();
                return false;
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of a complex number.
        /// </summary>
        /// <returns>Complex number</returns>
        object ICloneable.Clone()
        {
            return new Complex32(this.Real, this.Imag);
        }
        /// <summary>
        /// Creates a copy of a complex number.
        /// </summary>
        /// <returns>Complex number</returns>
        public Complex32 Clone()
        {
            return new Complex32(this.Real, this.Imag);
        }
        #endregion

        #region Static methods

        /// <summary>
        /// Creates a complex number from polar coordinates.
        /// </summary>
        /// <param name="magnitude">Magnitude (radius)</param>
        /// <param name="phase">Phase (angle in radians)</param>
        public static Complex32 FromPolarCoordinates(float magnitude, float phase)
        {
            return new Complex32(
                magnitude * (float)Math.Cos(phase),
                magnitude * (float)Math.Sin(phase)
            );
        }

        /// <summary>
        /// Returns a value that indicates whether the specified value is not a number.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Boolean</returns>
        public static bool IsNaN(Complex32 z) => float.IsNaN(z.Real) || float.IsNaN(z.Imag);

        /// <summary>
        /// Returns a value indicating whether the specified number evaluates to negative or positive infinity.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Boolean</returns>
        public static bool IsInfinity(Complex32 z) => float.IsInfinity(z.Real) || float.IsInfinity(z.Imag);

        #endregion

        #region Static properties

        /// <summary>
        /// Returns the imaginary one.
        /// </summary>
        public static Complex32 I
        {
            get
            {
                return new Complex32(0, 1);
            }
        }
        /// <summary>
        /// Returns the real one.
        /// </summary>
        public static Complex32 One
        {
            get
            {
                return new Complex32(1, 0);
            }
        }
        /// <summary>
        /// Returns the complex zero.
        /// </summary>
        public static Complex32 Zero
        {
            get
            {
                return new Complex32(0, 0);
            }
        }
        /// <summary>
        /// Returns the complex conjugate number.
        /// </summary>
        public Complex32 Conjugate
        {
            get
            {
                return new Complex32(this.Real, -this.Imag);
            }
        }
        /// <summary>
        /// Returns the not number.
        /// </summary>
        public static Complex32 NaN
        {
            get
            {
                return new Complex32(float.NaN, float.NaN);
            }
        }

        #endregion
    }
}
