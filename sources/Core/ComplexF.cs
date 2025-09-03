using System;
using System.Numerics;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a complex number.
    /// </summary>
    [Serializable]
    public struct ComplexF : ICloneable
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
        public ComplexF(float re, float im)
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
            return (obj is ComplexF) ? (this == (ComplexF)obj) : false;
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
        public static bool operator ==(ComplexF a, ComplexF b)
        {
            return ((a.Real == b.Real) && (a.Imag == b.Imag));
        }
        /// <summary>
        /// Checks if two complex numbers are not equal.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(ComplexF a, ComplexF b)
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
        public static ComplexF operator +(ComplexF a, ComplexF b)
        {
            return new ComplexF(a.Real + b.Real, a.Imag + b.Imag);
        }
        /// <summary>
        /// The sum of a complex number and a real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Value</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator +(ComplexF a, float b)
        {
            return new ComplexF(a.Real + b, a.Imag);
        }
        /// <summary>
        /// The sum of a complex number and a real number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator +(float a, ComplexF b)
        {
            return new ComplexF(b.Real + a, b.Imag);
        }


        /// <summary>
        /// The difference of two complex numbers.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator -(ComplexF a, ComplexF b)
        {
            return new ComplexF(a.Real - b.Real, a.Imag - b.Imag);
        }
        /// <summary>
        /// The difference between a complex number and a real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Value</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator -(ComplexF a, float b)
        {
            return new ComplexF(a.Real - b, a.Imag);
        }
        /// <summary>
        /// The difference between a complex number and a real number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator -(float a, ComplexF b)
        {
            return new ComplexF(a - b.Real, -b.Imag);
        }
        /// <summary>
        /// Inverts complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator -(ComplexF a)
        {
            return new ComplexF(-a.Real, -a.Imag);
        }


        /// <summary>
        /// Multiplies one complex number by another.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator *(ComplexF a, ComplexF b)
        {
            float aRe = a.Real, aIm = a.Imag;
            float bRe = b.Real, bIm = b.Imag;

            return new ComplexF(aRe * bRe - aIm * bIm, aRe * bIm + aIm * bRe);
        }
        /// <summary>
        /// Multiplies real number by complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Value</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator *(float a, ComplexF b)
        {
            return new ComplexF(b.Real * a, b.Imag * a);
        }
        /// <summary>
        /// Multiplies complex number by real number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator *(ComplexF a, float b)
        {
            return new ComplexF(a.Real * b, a.Imag * b);
        }


        /// <summary>
        /// Divides one complex number by another.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator /(ComplexF a, ComplexF b)
        {
            float aRe = a.Real, aIm = a.Imag;
            float bRe = b.Real, bIm = b.Imag;
            float abs = bRe * bRe + bIm * bIm;
            float inv = 1 / abs;

            return new ComplexF((aRe * bRe + aIm * bIm) * inv, (aIm * bRe - aRe * bIm) * inv);
        }
        /// <summary>
        /// Divides complex number by real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Value</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator /(ComplexF a, float b)
        {
            return new ComplexF(a.Real / b, a.Imag / b);
        }
        /// <summary>
        /// Divides real number by complex number.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static ComplexF operator /(float a, ComplexF b)
        {
            // (a + 0i) / (bRe + i*bIm) = a*(bRe - i*bIm) / (bRe^2 + bIm^2)
            float bRe = b.Real;
            float bIm = b.Imag;
            float abs = bRe * bRe + bIm * bIm;

            return new ComplexF(a * bRe / abs, -a * bIm / abs);
        }
        #endregion

        #region Conversion operators
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(Complex value)
        {
            return new ComplexF((float)value.Real, (float)value.Imaginary);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(ComplexF value)
        {
            return new Complex(value.Real, value.Imag);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(double value)
        {
            return new ComplexF((float)value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(float value)
        {
            return new ComplexF(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(long value)
        {
            return new ComplexF(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(ulong value)
        {
            return new ComplexF(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(short value)
        {
            return new ComplexF(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(ushort value)
        {
            return new ComplexF(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(int value)
        {
            return new ComplexF(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(uint value)
        {
            return new ComplexF(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(byte value)
        {
            return new ComplexF(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(sbyte value)
        {
            return new ComplexF(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator ComplexF(decimal value)
        {
            return new ComplexF((float)value, 0);
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
        public static ComplexF Parse(string s)
        {
            return StringOptions.Compar(s);
        }
        /// <summary>
        /// Tries to parse the string to complex number.
        /// </summary>
        /// <param name="complex">Input string</param>
        /// <param name="result">Complex number</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string complex, out ComplexF result)
        {
            try
            {
                result = ComplexF.Parse(complex);
                return true;
            }
            catch (FormatException)
            {
                result = new ComplexF();
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
            return new ComplexF(this.Real, this.Imag);
        }
        /// <summary>
        /// Creates a copy of a complex number.
        /// </summary>
        /// <returns>Complex number</returns>
        public ComplexF Clone()
        {
            return new ComplexF(this.Real, this.Imag);
        }
        #endregion

        #region Static methods

        /// <summary>
        /// Creates a complex number from polar coordinates.
        /// </summary>
        /// <param name="magnitude">Magnitude (radius)</param>
        /// <param name="phase">Phase (angle in radians)</param>
        public static ComplexF FromPolarCoordinates(float magnitude, float phase)
        {
            return new ComplexF(
                magnitude * (float)Math.Cos(phase),
                magnitude * (float)Math.Sin(phase)
            );
        }

        /// <summary>
        /// Returns a value that indicates whether the specified value is not a number.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Boolean</returns>
        public static bool IsNaN(ComplexF z) => float.IsNaN(z.Real) || float.IsNaN(z.Imag);

        /// <summary>
        /// Returns a value indicating whether the specified number evaluates to negative or positive infinity.
        /// </summary>
        /// <param name="z">Value</param>
        /// <returns>Boolean</returns>
        public static bool IsInfinity(ComplexF z) => float.IsInfinity(z.Real) || float.IsInfinity(z.Imag);

        #endregion

        #region Static properties

        /// <summary>
        /// Returns the imaginary one.
        /// </summary>
        public static ComplexF I
        {
            get
            {
                return new ComplexF(0, 1);
            }
        }
        /// <summary>
        /// Returns the real one.
        /// </summary>
        public static ComplexF One
        {
            get
            {
                return new ComplexF(1, 0);
            }
        }
        /// <summary>
        /// Returns the complex zero.
        /// </summary>
        public static ComplexF Zero
        {
            get
            {
                return new ComplexF(0, 0);
            }
        }
        /// <summary>
        /// Returns the complex conjugate number.
        /// </summary>
        public ComplexF Conjugate
        {
            get
            {
                return new ComplexF(this.Real, -this.Imag);
            }
        }
        /// <summary>
        /// Returns the not number.
        /// </summary>
        public static ComplexF NaN
        {
            get
            {
                return new ComplexF(float.NaN, float.NaN);
            }
        }

        #endregion
    }
}
