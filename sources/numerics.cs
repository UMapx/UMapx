// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Runtime.Serialization;
using System.Text.RegularExpressions;

namespace UMapx.Core
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                                 UMAPX.CORE
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Complex
    /// <summary>
    /// Defines a complex number.
    /// </summary>
    public struct Complex : ICloneable, ISerializable
    {
        #region Private data
        /// <summary>
        /// The real part of the complex number.
        /// </summary>
        public double Real;
        /// <summary>
        /// The imaginary part of a complex number.
        /// </summary>
        public double Imag;
        #endregion

        #region Structure components
        /// <summary>
        /// Initializes the complex number.
        /// </summary>
        /// <param name="re">Real part of the complex number</param>
        /// <param name="im">Imaginary part of a complex number</param>
        public Complex(double re, double im)
        {
            this.Real = re;
            this.Imag = im;
        }
        /// <summary>
        /// Gets the value of the module.
        /// </summary>
        public double Abs
        {
            get
            {
                return Math.Sqrt(Real * Real + Imag * Imag);
            }
        }
        /// <summary>
        /// Gets the value of the phase.
        /// </summary>
        public double Angle
        {
            get
            {
                return Math.Atan2(Imag, Real);
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
            return StringOptions.Disp(new double[] { this.Real, this.Imag }, format, StringOptions.C);
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
        public static Complex operator +(Complex a, double b)
        {
            return new Complex(a.Real + b, a.Imag);
        }
        /// <summary>
        /// The sum of a complex number and a real number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator +(double a, Complex b)
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
        public static Complex operator -(Complex a, double b)
        {
            return new Complex(a.Real - b, a.Imag);
        }
        /// <summary>
        /// The difference between a complex number and a real number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator -(double a, Complex b)
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
            double aRe = a.Real, aIm = a.Imag;
            double bRe = b.Real, bIm = b.Imag;

            return new Complex(aRe * bRe - aIm * bIm, aRe * bIm + aIm * bRe);
        }
        /// <summary>
        /// Multiplies real number by complex number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Number</param>
        /// <returns>Complex number</returns>
        public static Complex operator *(double a, Complex b)
        {
            return new Complex(b.Real * a, b.Imag * a);
        }
        /// <summary>
        /// Multiplies complex number by real number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator *(Complex a, double b)
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
            double aRe = a.Real, aIm = a.Imag;
            double bRe = b.Real, bIm = b.Imag;
            double abs = bRe * bRe + bIm * bIm;
            double inv = 1 / abs;

            return new Complex((aRe * bRe + aIm * bIm) * inv, (aIm * bRe - aRe * bIm) * inv);
        }
        /// <summary>
        /// Divides complex number by real number.
        /// </summary>
        /// <param name="a">Complex number</param>
        /// <param name="b">Number</param>
        /// <returns>Complex number</returns>
        public static Complex operator /(Complex a, double b)
        {
            return new Complex(a.Real / b, a.Imag / b);
        }
        /// <summary>
        /// Divides real number by complex number.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="b">Complex number</param>
        /// <returns>Complex number</returns>
        public static Complex operator /(double a, Complex b)
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
            return new Complex(value, 0);
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
            return new Complex((double)value, 0);
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
        public static bool TryParse(string complex, ref Complex result)
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

        #region Serialization members
        /// <summary>
        /// Gets information about an object.
        /// </summary>
        /// <param name="info">Data needed for object serialization and deserialization</param>
        /// <param name="context">Source and destination of a given stream</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("Real", this.Real);
            info.AddValue("Imaginary", this.Imag);
        }
        #endregion
    }
    #endregion

    #region Quaternion
    /// <summary>
    /// Defines a quaternion.
    /// <remarks>
    /// A quaternion is a system of hypercomplex numbers that forms a four-dimensional vector space over a field of real numbers.
    /// </remarks>
    /// </summary>
    public struct Quaternion : ICloneable, ISerializable
    {
        #region Public data
        /// <summary>
        /// X coordinate.
        /// </summary>
        public double X;
        /// <summary>
        /// Y coordinate.
        /// </summary>
        public double Y;
        /// <summary>
        /// Z coordinate.
        /// </summary>
        public double Z;
        /// <summary>
        /// W coordinate.
        /// </summary>
        public double W;
        #endregion

        #region Quaternion components
        /// <summary>
        /// Creates a quaternion based on the given coordinates.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        /// <param name="z">Coordinate Z</param>
        /// <param name="w">Coordinate W</param>
        public Quaternion(double x, double y, double z, double w)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
            this.W = w;
        }
        /// <summary>
        /// Gets a quaternion that represents a lack of rotation.
        /// </summary>
        public static Quaternion Identity
        {

            get
            {
                return new Quaternion(0, 0, 0, 1);
            }
        }
        /// <summary>
        /// Gets a value indicating whether the current instance is a single Quaternion.
        /// </summary>
        public bool IsIdentity
        {
            get
            {
                return this.X == 0.0 && this.Y == 0.0 && this.Z == 0.0 && this.W == 1.0;
            }
        }
        /// <summary>
        /// Returns the value of the quaternion module.
        /// </summary>
        public double Abs
        {
            get
            {
                return Math.Sqrt(this.X * this.X + this.Y * this.Y + this.Z * this.Z + this.W * this.W);
            }
        }
        /// <summary>
        /// Calculates the quaternion modulus squared.
        /// </summary>
        public double SquaredAbs
        {
            get
            {
                return this.X * this.X + this.Y * this.Y + this.Z * this.Z + this.W * this.W;
            }
        }
        #endregion

        #region Operations
        /// <summary>
        /// Divides each coordinate of the specified quaternion by its length.
        /// </summary>
        public Quaternion Normalize
        {
            get
            {
                double norm = 1.0 / this.Abs;
                return new Quaternion(this.X * norm, this.Y * norm, this.Z * norm, this.W * norm);
            }
        }
        /// <summary>
        /// Returns the conjugate object of the specified quaternion.
        /// </summary>
        public Quaternion Conjugate
        {
            get
            {
                return new Quaternion(-this.X, -this.Y, -this.Z, this.W);
            }
        }
        /// <summary>
        /// Returns the inverse object of the quaternion.
        /// </summary>
        public Quaternion Inverse
        {
            get
            {
                double norm = 1.0 / this.SquaredAbs;
                return new Quaternion(-this.X * norm, -this.Y * norm, -this.Z * norm, this.W * norm);
            }
        }
        /// <summary>
        /// Creates a new quaternion based on a given value of nutation, precession, and proper rotation.
        /// </summary>
        /// <param name="yaw">The nutation angle around the Y axis in radians</param>
        /// <param name="pitch">The precession angle around the X axis in radians</param>
        /// <param name="roll">The angle of rotation around the Z axis in radians</param>
        /// <returns>Quaternion</returns>
        public static Quaternion FromYPR(double yaw, double pitch, double roll)
        {
            double a = roll * 0.5;
            double b = Math.Sin(a);
            double c = Math.Cos(a);
            double d = pitch * 0.5;
            double e = Math.Sin(d);
            double f = Math.Cos(d);
            double g = yaw * 0.5;
            double h = Math.Sin(g);
            double i = Math.Cos(g);

            return new Quaternion(
                i * e * c + h * f * b,
                h * f * c - i * e * b,
                i * f * b - h * e * c,
                i * f * c + h * e * b);
        }
        /// <summary>
        /// Computes the scalar product of two quaternion.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Quaternion</returns>
        public static double Dot(Quaternion a, Quaternion b)
        {
            return a.X * b.X + a.Y * b.Y + a.Z * b.Z + a.W * b.W;
        }
        /// <summary>
        /// Performs interpolation between two quaternions using spherical linear interpolation.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <param name="amount">Relative weight of the second quaternion in interpolation</param>
        /// <returns>Quaternion</returns>
        public static Quaternion Slerp(Quaternion a, Quaternion b, double amount)
        {
            double d, e, dot = Quaternion.Dot(a, b);
            bool flag = false;

            if (dot < 0.0)
            {
                flag = true;
                dot = -dot;
            }
            if (dot > 0.999999)
            {
                d = 1.0 - amount;
                e = (flag ? (-amount) : amount);
            }
            else
            {
                double f = Math.Acos(dot);
                double g = (1.0 / Math.Sin(f));
                d = Math.Sin(((1.0 - amount) * f)) * g;
                e = (flag ? ((-Math.Sin((amount * f))) * g) : (Math.Sin((amount * f)) * g));
            }

            return new Quaternion(
                d * a.X + e * b.X,
                d * a.Y + e * b.Y,
                d * a.Z + e * b.Z,
                d * a.W + e * b.W);
        }
        /// <summary>
        /// Performs linear interpolation between two quaternions based on a value indicating the weighting of the second quaternion.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <param name="amount">Relative weight of the second quaternion in interpolation</param>
        /// <returns>Quaternion</returns>
        public static Quaternion Lerp(Quaternion a, Quaternion b, double amount)
        {
            double f = 1.0 - amount;
            Quaternion quaternion3 = default(Quaternion);
            double dot = Dot(a, b);

            if (dot >= 0.0)
            {
                quaternion3.X = f * a.X + amount * b.X;
                quaternion3.Y = f * a.Y + amount * b.Y;
                quaternion3.Z = f * a.Z + amount * b.Z;
                quaternion3.W = f * a.W + amount * b.W;
            }
            else
            {
                quaternion3.X = f * a.X - amount * b.X;
                quaternion3.Y = f * a.Y - amount * b.Y;
                quaternion3.Z = f * a.Z - amount * b.Z;
                quaternion3.W = f * a.W - amount * b.W;
            }

            double norm = 1.0 / quaternion3.Abs;
            quaternion3.X *= norm;
            quaternion3.Y *= norm;
            quaternion3.Z *= norm;
            quaternion3.W *= norm;
            return quaternion3;
        }
        /// <summary>
        /// Concatenates two quaternions.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Quaternion</returns>
        public static Quaternion Concatenate(Quaternion a, Quaternion b)
        {
            double x = b.X, y = b.Y, z = b.Z, w = b.W;
            double x2 = a.X, y2 = a.Y, z2 = a.Z, w2 = a.W;

            double e = y * z2 - z * y2;
            double f = z * x2 - x * z2;
            double c = x * y2 - y * x2;
            double d = x * x2 + y * y2 + z * z2;

            return new Quaternion(
                x * w2 + x2 * w + e,
                y * w2 + y2 * w + f,
                z * w2 + z2 * w + c,
                w * w2 - d);
        }
        #endregion

        #region Operators
        /// <summary>
        /// Reverses the sign of each quaternion coordinate.
        /// </summary>
        /// <param name="q">Quaternion</param>
        /// <returns>Quaternion</returns>
        public static Quaternion operator -(Quaternion q)
        {
            return new Quaternion(-q.X, -q.Y, -q.Z, -q.W);
        }
        /// <summary>
        /// Adds each element in one quaternion with the corresponding element in the second quaternion.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Quaternion</returns>
        public static Quaternion operator +(Quaternion a, Quaternion b)
        {
            return new Quaternion(
                a.X + b.X,
                a.Y + b.Y,
                a.Z + b.Z,
                a.W + b.W);
        }
        /// <summary>
        /// Subtracts each element in the second quaternion from the corresponding element in the first quaternion.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Quaternion</returns>
        public static Quaternion operator -(Quaternion a, Quaternion b)
        {
            return new Quaternion(
                a.X - b.X,
                a.Y - b.Y,
                a.Z - b.Z,
                a.W - b.W);
        }
        /// <summary>
        /// Returns the quaternion resulting from the multiplication of two quaternions.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Quaternion</returns>
        public static Quaternion operator *(Quaternion a, Quaternion b)
        {
            double x = a.X;
            double y = a.Y;
            double z = a.Z;
            double w = a.W;
            double x2 = b.X;
            double y2 = b.Y;
            double z2 = b.Z;
            double w2 = b.W;
            double d = y * z2 - z * y2;
            double e = z * x2 - x * z2;
            double g = x * y2 - y * x2;
            double h = x * x2 + y * y2 + z * z2;

            return new Quaternion(
                 x * w2 + x2 * w + d,
                 y * w2 + y2 * w + e,
                 z * w2 + z2 * w + g,
                 w * w2 - h
                );
        }
        /// <summary>
        /// Returns the quaternion obtained by scaling all the coordinates of the specified quaternion by a scalar factor.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Factor</param>
        /// <returns>Quaternion</returns>
        public static Quaternion operator *(Quaternion a, double b)
        {
            return new Quaternion(
                a.X * b,
                a.Y * b,
                a.Z * b,
                a.W * b);
        }
        /// <summary>
        /// Divides one quaternion into a second quaternion.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Quaternion</returns>
        public static Quaternion operator /(Quaternion a, Quaternion b)
        {
            double x = a.X;
            double y = a.Y;
            double z = a.Z;
            double w = a.W;
            double d = 1.0 / b.SquaredAbs;
            double e = -b.X * d;
            double f = -b.Y * d;
            double g = -b.Z * d;
            double i = b.W * d;
            double j = y * g - z * f;
            double k = z * e - x * g;
            double l = x * f - y * e;
            double m = x * e + y * f + z * g;

            return new Quaternion(
                x * i + e * w + j,
                y * i + f * w + k,
                z * i + g * w + l,
                w * i - m);
        }
        /// <summary>
        /// Returns the quaternion obtained by scaling all the coordinates of the specified quaternion by a scalar factor.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Factor</param>
        /// <returns>Quaternion</returns>
        public static Quaternion operator /(Quaternion a, double b)
        {
            return new Quaternion(
                a.X / b,
                a.Y / b,
                a.Z / b,
                a.W / b);
        }
        #endregion

        #region Bools & overrides
        /// <summary>
        /// Checks if two quaternions are equal.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(Quaternion a, Quaternion b)
        {
            return a.X == b.X && a.Y == b.Y && a.Z == b.Z && a.W == b.W;
        }
        /// <summary>
        /// Checks if two quaternions are not equal.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(Quaternion a, Quaternion b)
        {
            return !(a == b);
        }
        /// <summary>
        /// Gets a value indicating whether this instance is equal to the specified value of type quaternion.
        /// </summary>
        /// <param name="obj">Object</param>
        /// <returns>Boolean</returns>
        public override bool Equals(object obj)
        {
            return (obj is Quaternion) ? (this == (Quaternion)obj) : false;
        }
        /// <summary>
        /// Converts quaternion to its corresponding string representation.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return this.ToString("G3");
        }
        /// <summary>
        /// Converts quaternion to its corresponding string representation.
        /// </summary>
        /// <param name="format">Format string</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public string ToString(string format)
        {
            return StringOptions.Disp(new double[] { this.X, this.Y, this.Z, this.W }, format, StringOptions.Q);
        }
        /// <summary>
        /// Returns the hash code for this object.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return this.X.GetHashCode() + this.Y.GetHashCode() + this.Z.GetHashCode() + this.W.GetHashCode();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of quaternion.
        /// </summary>
        /// <returns>Quaternion</returns>
        object ICloneable.Clone()
        {
            return new Quaternion(this.X, this.Y, this.Z, this.W);
        }
        /// <summary>
        /// Creates a copy of quaternion.
        /// </summary>
        /// <returns>Quaternion</returns>
        public Quaternion Clone()
        {
            return new Quaternion(this.X, this.Y, this.Z, this.W);
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Gets information about an object.
        /// </summary>
        /// <param name="info">Data needed for object serialization and deserialization</param>
        /// <param name="context">Source and destination of a given stream</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("X", this.X);
            info.AddValue("Y", this.Y);
            info.AddValue("Z", this.Z);
            info.AddValue("W", this.W);
        }
        #endregion

        #region Parsing
        /// <summary>
        /// Parses the string to quaternion.
        /// </summary>
        /// <remarks>
        /// Example: "[1, -2; 3.2, -.13]";
        /// The dimension of the vector must be 4.
        /// </remarks>
        /// <param name="s">Input string</param>
        /// <returns>Quaternion</returns>
        public static Quaternion Parse(string s)
        {
            string[] cols = StringOptions.Matpar(s);
            string[] nums = cols[0].Split('|');

            if (cols.Length > 1 || nums.Length != 4)
                throw new Exception("The input string was in the wrong format");

            return new Quaternion(double.Parse(nums[0]),
                                  double.Parse(nums[1]),
                                  double.Parse(nums[2]),
                                  double.Parse(nums[3]));
        }
        /// <summary>
        /// Tries to parse the string into Quaternion.
        /// </summary>
        /// <param name="quaternion">Input string</param>
        /// <param name="result">Quaternion</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string quaternion, ref Quaternion result)
        {
            try
            {
                result = Quaternion.Parse(quaternion);
                return true;
            }
            catch (FormatException)
            {
                result = new Quaternion();
                return false;
            }
        }
        #endregion
    }
    #endregion

    #region Internal class
    /// <summary>
    /// Defines a class of string operations.
    /// </summary>
    internal class StringOptions
    {
        #region String voids
        /// <summary>
        /// Complex number format.
        /// </summary>
        public static string[] C
        {
            get
            {
                return new string[] { "", "i" };
            }
        }
        /// <summary>
        /// Quaternion format.
        /// </summary>
        public static string[] Q
        {
            get
            {
                return new string[] { "i", "j", "k", "" };
            }
        }
        /// <summary>
        /// The function of converting an array of numbers to a string.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="format">Format string</param>
        /// <param name="symbol">String array</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Disp(double[] v, string format, string[] symbol)
        {
            int length = v.Length, i;
            int start = -1;

            for (i = 0; i < length; i++)
            {
                if (v[i] != 0)
                {
                    start = i;
                    break;
                }
            }

            if (start != -1)
            {
                string result = Disp(v[start], format, true, symbol[start]);

                for (i = start + 1; i < length; i++)
                {
                    result += Disp(v[i], format, false, symbol[i]);
                }
                return result;
            }
            return "0";
        }
        /// <summary>
        /// The function of converting number to a string
        /// </summary>
        /// <param name="v">Number</param>
        /// <param name="format">Format string</param>
        /// <param name="s">First in a row or not</param>
        /// <param name="symbol">Symbol</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Disp(double v, string format, bool s, string symbol)
        {
            if (v == 0)
            {
                return "";
            }
            else if (v < 0)
            {
                return (s) ? "-" + (-v).ToString(format) + symbol : " - " + (-v).ToString(format) + symbol;
            }
            return (s) ? v.ToString(format) + symbol : " + " + v.ToString(format) + symbol;
        }
        /// <summary>
        /// Defines a general method for casting the original row to the matrix form.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <returns>String array</returns>
        public static string[] Matpar(string s)
        {
            // example: s = "[1,2,3,4]".
            // Regex options:
            Regex regex = new Regex(@"\[(?<matrice>.*)]", RegexOptions.None);
            Match match = regex.Match(s);

            // success?
            if (match.Success)
            {
                // get new string:
                return Regex.Split(match.Result("${matrice}").Replace(",", "|"), ";");
            }
            throw new Exception("The input string was in the wrong format");
        }
        /// <summary>
        /// Translates the original string to complex number.
        /// <remarks>
        /// Example: "1 + 2i", "0.321 + 11i", ".1i".
        /// </remarks>
        /// </summary>
        /// <param name="s">Input string</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static Complex Compar(string s)
        {
            string u = s.Replace(" ", "");
            int i, k = 0;            
            int length = u.Length;     
            string re = "" + u[0];       
            string im = "";            

            if (!u.Contains("i"))
            {
                for (i = 1; i < length; i++)
                {
                    if (u[i] != '+' && u[i] != '-')
                    {
                        re += u[i]; k++;
                    }
                    else break;
                }

                return new Complex(double.Parse(re), 0);
            }
            else
            {
                if (u == "i") return new Complex(0, 1);

                for (i = 1; i < length; i++)
                {
                    if (u[i] != '+' && u[i] != '-')
                    {
                        re += u[i]; k++;
                    }
                    else break;
                }

                if (k != length - 1)
                {
                    int k1 = k + 1, k2 = k + 2; im += u[k1];

                    if (u[k2] == 'i') return new Complex(double.Parse(re), double.Parse(im + '1'));

                    for (i = k2; i < length; i++)
                    {
                        if (u[i] != 'i')
                        {
                            im += u[i];
                        }
                        else break;
                    }
                    return new Complex(double.Parse(re), double.Parse(im));
                }
                else
                {
                    return new Complex(0, double.Parse(re.Replace("i", "")));
                }
            }
        }
        #endregion
    }
    #endregion
}
