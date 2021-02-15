using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a complex number.
    /// </summary>
    [Serializable]
    public struct Complex : ICloneable
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
        public Complex(float re, float im)
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
        /// Gets the value of the phase.
        /// </summary>
        public float Angle
        {
            get
            {
                return (float)Math.Atan2(Imag, Real);
            }
        }
        /// <summary>
        /// Returns the imaginary one.
        /// </summary>
        public static Complex I
        {
            get
            {
                return new Complex(0, 1);
            }
        }
        /// <summary>
        /// Returns the real one.
        /// </summary>
        public static Complex One
        {
            get
            {
                return new Complex(1, 0);
            }
        }
        /// <summary>
        /// Returns the complex zero.
        /// </summary>
        public static Complex Zero
        {
            get
            {
                return new Complex(0, 0);
            }
        }
        /// <summary>
        /// Returns the complex conjugate number.
        /// </summary>
        public Complex Conjugate
        {
            get
            {
                return new Complex(this.Real, -this.Imag);
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
            return (obj is Complex) ? (this == (Complex)obj) : false;
        }
        /// <summary>
        /// Converts complex number to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return this.ToString("G3");
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
        public static bool operator ==(Complex a, Complex b)
        {
            return ((a.Real == b.Real) && (a.Imag == b.Imag));
        }
        /// <summary>
        /// Checks if two complex numbers are not equal.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(Complex a, Complex b)
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
        public static Complex operator +(Complex a, Complex b)
        {
            return new Complex(a.Real + b.Real, a.Imag + b.Imag);
        }
        /// <summary>
        /// The sum of a complex number and a real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Number</param>
        /// <returns>Complex number</returns>
        public static Complex operator +(Complex a, float b)
        {
            return new Complex(a.Real + b, a.Imag);
        }
        /// <summary>
        /// The sum of a complex number and a real number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator +(float a, Complex b)
        {
            return new Complex(b.Real + a, b.Imag);
        }


        /// <summary>
        /// The difference of two complex numbers.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator -(Complex a, Complex b)
        {
            return new Complex(a.Real - b.Real, a.Imag - b.Imag);
        }
        /// <summary>
        /// The difference between a complex number and a real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Number</param>
        /// <returns>Complex number</returns>
        public static Complex operator -(Complex a, float b)
        {
            return new Complex(a.Real - b, a.Imag);
        }
        /// <summary>
        /// The difference between a complex number and a real number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator -(float a, Complex b)
        {
            return new Complex(a - b.Real, b.Imag);
        }
        /// <summary>
        /// Inverts complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator -(Complex a)
        {
            return new Complex(-a.Real, -a.Imag);
        }


        /// <summary>
        /// Multiplies one complex number by another.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator *(Complex a, Complex b)
        {
            float aRe = a.Real, aIm = a.Imag;
            float bRe = b.Real, bIm = b.Imag;

            return new Complex(aRe * bRe - aIm * bIm, aRe * bIm + aIm * bRe);
        }
        /// <summary>
        /// Multiplies real number by complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Number</param>
        /// <returns>Complex number</returns>
        public static Complex operator *(float a, Complex b)
        {
            return new Complex(b.Real * a, b.Imag * a);
        }
        /// <summary>
        /// Multiplies complex number by real number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator *(Complex a, float b)
        {
            return new Complex(a.Real * b, a.Imag * b);
        }


        /// <summary>
        /// Divides one complex number by another.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator /(Complex a, Complex b)
        {
            float aRe = a.Real, aIm = a.Imag;
            float bRe = b.Real, bIm = b.Imag;
            float abs = bRe * bRe + bIm * bIm;
            float inv = 1 / abs;

            return new Complex((aRe * bRe + aIm * bIm) * inv, (aIm * bRe - aRe * bIm) * inv);
        }
        /// <summary>
        /// Divides complex number by real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Number</param>
        /// <returns>Complex number</returns>
        public static Complex operator /(Complex a, float b)
        {
            return new Complex(a.Real / b, a.Imag / b);
        }
        /// <summary>
        /// Divides real number by complex number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator /(float a, Complex b)
        {
            if (b.Imag == 0)
            {
                return new Complex(a / b.Real, 0);
            }
            else if (b.Real == 0)
            {
                return new Complex(0, a / b.Imag);
            }
            return new Complex(a / b.Real, a / b.Imag);
        }
        #endregion

        #region Conversion operators
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(double value)
        {
            return new Complex((float)value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(float value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(long value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(ulong value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(short value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(ushort value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(int value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(uint value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(byte value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(sbyte value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Defines an explicit conversion of a number to complex number.
        /// </summary>
        /// <param name="value">Value to be converted to complex number</param>
        /// <returns>Complex number</returns>
        public static implicit operator Complex(decimal value)
        {
            return new Complex((float)value, 0);
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
        public static Complex Parse(string s)
        {
            return StringOptions.Compar(s);
        }
        /// <summary>
        /// Tries to parse the string to complex number.
        /// </summary>
        /// <param name="complex">Input string</param>
        /// <param name="result">Complex number</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string complex, out Complex result)
        {
            try
            {
                result = Complex.Parse(complex);
                return true;
            }
            catch (FormatException)
            {
                result = new Complex();
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
            return new Complex(this.Real, this.Imag);
        }
        /// <summary>
        /// Creates a copy of a complex number.
        /// </summary>
        /// <returns>Complex number</returns>
        public Complex Clone()
        {
            return new Complex(this.Real, this.Imag);
        }
        #endregion
    }
}
