using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a quaternion.
    /// </summary>
    /// <remarks>
    /// A quaternion is a system of hypercomplex numbers that forms a four-dimensional vector space over a field of real numbers.
    /// </remarks>
    [Serializable]
    public struct Quaternion32 : ICloneable
    {
        #region Public data
        /// <summary>
        /// X coordinate.
        /// </summary>
        public float X;
        /// <summary>
        /// Y coordinate.
        /// </summary>
        public float Y;
        /// <summary>
        /// Z coordinate.
        /// </summary>
        public float Z;
        /// <summary>
        /// W coordinate.
        /// </summary>
        public float W;
        #endregion

        #region Quaternion components
        /// <summary>
        /// Creates a quaternion based on the given coordinates.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        /// <param name="z">Coordinate Z</param>
        /// <param name="w">Coordinate W</param>
        public Quaternion32(float x, float y, float z, float w)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
            this.W = w;
        }
        /// <summary>
        /// Gets a quaternion that represents a lack of rotation.
        /// </summary>
        public static Quaternion32 Identity
        {

            get
            {
                return new Quaternion32(0, 0, 0, 1);
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
        /// Returns the value of the quaternion modulus.
        /// </summary>
        public float Abs
        {
            get
            {
                return (float)Math.Sqrt(this.X * this.X + this.Y * this.Y + this.Z * this.Z + this.W * this.W);
            }
        }
        /// <summary>
        /// Calculates the quaternion modulus squared.
        /// </summary>
        public float SquaredAbs
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
        public Quaternion32 Normalize
        {
            get
            {
                float norm = 1.0f / this.Abs;
                return new Quaternion32(this.X * norm, this.Y * norm, this.Z * norm, this.W * norm);
            }
        }
        /// <summary>
        /// Returns the conjugate object of the specified quaternion.
        /// </summary>
        public Quaternion32 Conjugate
        {
            get
            {
                return new Quaternion32(-this.X, -this.Y, -this.Z, this.W);
            }
        }
        /// <summary>
        /// Returns the inverse object of the quaternion.
        /// </summary>
        public Quaternion32 Inverse
        {
            get
            {
                float norm = 1.0f / this.SquaredAbs;
                return new Quaternion32(-this.X * norm, -this.Y * norm, -this.Z * norm, this.W * norm);
            }
        }
        /// <summary>
        /// Creates a new quaternion based on a given value of nutation, precession, and proper rotation.
        /// </summary>
        /// <param name="yaw">The nutation angle around the Y axis in radians</param>
        /// <param name="pitch">The precession angle around the X axis in radians</param>
        /// <param name="roll">The angle of rotation around the Z axis in radians</param>
        /// <returns>Quaternion</returns>
        public static Quaternion32 FromYPR(float yaw, float pitch, float roll)
        {
            float a = roll * 0.5f;
            float b = (float)Math.Sin(a);
            float c = (float)Math.Cos(a);
            float d = pitch * 0.5f;
            float e = (float)Math.Sin(d);
            float f = (float)Math.Cos(d);
            float g = yaw * 0.5f;
            float h = (float)Math.Sin(g);
            float i = (float)Math.Cos(g);

            return new Quaternion32(
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
        public static float Dot(Quaternion32 a, Quaternion32 b)
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
        public static Quaternion32 Slerp(Quaternion32 a, Quaternion32 b, float amount)
        {
            float d, e, dot = Quaternion32.Dot(a, b);
            bool flag = false;

            if (dot < 0.0)
            {
                flag = true;
                dot = -dot;
            }
            if (dot > 0.999999)
            {
                d = 1.0f - amount;
                e = (flag ? (-amount) : amount);
            }
            else
            {
                float f = (float)Math.Acos(dot);
                float g = (float)(1.0 / Math.Sin(f));
                d = (float)Math.Sin(((1.0 - amount) * f)) * g;
                e = (float)(flag ? ((-Math.Sin((amount * f))) * g) : (Math.Sin((amount * f)) * g));
            }

            return new Quaternion32(
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
        public static Quaternion32 Lerp(Quaternion32 a, Quaternion32 b, float amount)
        {
            float f = 1.0f - amount;
            Quaternion32 quaternion3 = default(Quaternion32);
            float dot = Dot(a, b);

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

            float norm = 1.0f / quaternion3.Abs;
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
        public static Quaternion32 Concatenate(Quaternion32 a, Quaternion32 b)
        {
            float x = b.X, y = b.Y, z = b.Z, w = b.W;
            float x2 = a.X, y2 = a.Y, z2 = a.Z, w2 = a.W;

            float e = y * z2 - z * y2;
            float f = z * x2 - x * z2;
            float c = x * y2 - y * x2;
            float d = x * x2 + y * y2 + z * z2;

            return new Quaternion32(
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
        public static Quaternion32 operator -(Quaternion32 q)
        {
            return new Quaternion32(-q.X, -q.Y, -q.Z, -q.W);
        }
        /// <summary>
        /// Adds each element in one quaternion with the corresponding element in the second quaternion.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Quaternion</returns>
        public static Quaternion32 operator +(Quaternion32 a, Quaternion32 b)
        {
            return new Quaternion32(
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
        public static Quaternion32 operator -(Quaternion32 a, Quaternion32 b)
        {
            return new Quaternion32(
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
        public static Quaternion32 operator *(Quaternion32 a, Quaternion32 b)
        {
            float x = a.X;
            float y = a.Y;
            float z = a.Z;
            float w = a.W;
            float x2 = b.X;
            float y2 = b.Y;
            float z2 = b.Z;
            float w2 = b.W;
            float d = y * z2 - z * y2;
            float e = z * x2 - x * z2;
            float g = x * y2 - y * x2;
            float h = x * x2 + y * y2 + z * z2;

            return new Quaternion32(
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
        public static Quaternion32 operator *(Quaternion32 a, float b)
        {
            return new Quaternion32(
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
        public static Quaternion32 operator /(Quaternion32 a, Quaternion32 b)
        {
            float x = a.X;
            float y = a.Y;
            float z = a.Z;
            float w = a.W;
            float d = 1.0f / b.SquaredAbs;
            float e = -b.X * d;
            float f = -b.Y * d;
            float g = -b.Z * d;
            float i = b.W * d;
            float j = y * g - z * f;
            float k = z * e - x * g;
            float l = x * f - y * e;
            float m = x * e + y * f + z * g;

            return new Quaternion32(
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
        public static Quaternion32 operator /(Quaternion32 a, float b)
        {
            return new Quaternion32(
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
        public static bool operator ==(Quaternion32 a, Quaternion32 b)
        {
            return a.X == b.X && a.Y == b.Y && a.Z == b.Z && a.W == b.W;
        }
        /// <summary>
        /// Checks if two quaternions are not equal.
        /// </summary>
        /// <param name="a">Quaternion</param>
        /// <param name="b">Quaternion</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(Quaternion32 a, Quaternion32 b)
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
            return (obj is Quaternion32) ? (this == (Quaternion32)obj) : false;
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
            return StringOptions.Disp(new float[] { this.X, this.Y, this.Z, this.W }, format, StringOptions.Q);
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
            return new Quaternion32(this.X, this.Y, this.Z, this.W);
        }
        /// <summary>
        /// Creates a copy of quaternion.
        /// </summary>
        /// <returns>Quaternion</returns>
        public Quaternion32 Clone()
        {
            return new Quaternion32(this.X, this.Y, this.Z, this.W);
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
        public static Quaternion32 Parse(string s)
        {
            string[] cols = StringOptions.Matpar(s);
            string[] nums = cols[0].Split('|');

            if (cols.Length > 1 || nums.Length != 4)
                throw new ArgumentException("The input string was in the wrong format");

            return new Quaternion32(float.Parse(nums[0]),
                                  float.Parse(nums[1]),
                                  float.Parse(nums[2]),
                                  float.Parse(nums[3]));
        }
        /// <summary>
        /// Tries to parse the string into Quaternion.
        /// </summary>
        /// <param name="quaternion">Input string</param>
        /// <param name="result">Quaternion</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string quaternion, ref Quaternion32 result)
        {
            try
            {
                result = Quaternion32.Parse(quaternion);
                return true;
            }
            catch (FormatException)
            {
                result = new Quaternion32();
                return false;
            }
        }
        #endregion
    }
}
