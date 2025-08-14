using System;
using System.Drawing;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines an unknown color model.
    /// This color model can play the role of any color space.
    /// </summary>
	[Serializable]
    public struct Unknown : IColorSpace, ICloneable
    {
        #region Private data
        private float x;
        private float y;
        private float z;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure.
        /// </summary>
        /// <param name="x">Component X</param>
        /// <param name="y">Component Y</param>
        /// <param name="z">Component Z</param>
        public Unknown(float x, float y, float z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        /// <summary>
        /// Defines the component of the color model.
        /// </summary>
        public float X
        {
            get
            {
                return x;
            }
            set
            {
                x = value;
            }
        }
        /// <summary>
        /// Defines the component of the color model.
        /// </summary>
        public float Y
        {
            get
            {
                return y;
            }
            set
            {
                y = value;
            }
        }
        /// <summary>
        /// Defines the component of the color model.
        /// </summary>
        public float Z
        {
            get
            {
                return z;
            }
            set
            {
                z = value;
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">Unknown structure</param>
        /// <param name="item2">Unknown structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(Unknown item1, Unknown item2)
        {
            return (
                item1.X == item2.X
                && item1.Y == item2.Y
                && item1.Z == item2.Z
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">Unknown structure</param>
        /// <param name="item2">Unknown structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(Unknown item1, Unknown item2)
        {
            return !(item1 == item2);
        }
        #endregion

        #region Methods
        /// <summary>
        /// Defines whether the specified System.Object is equal to the current System.Object.
        /// </summary>
        /// <param name="obj">Element</param>
        /// <returns>Boolean</returns>
        public override bool Equals(Object obj)
        {
            if (obj == null || GetType() != obj.GetType()) return false;

            return (this == (Unknown)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return X.GetHashCode() ^ Y.GetHashCode() ^ Z.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return X.ToString() + Environment.NewLine + Y.ToString() + Environment.NewLine + Z.ToString();
        }
        #endregion

        #region Conversion operators
        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(AHSL value)
        {
            return new Unknown(value.Hue, value.Saturation, value.Lightness);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator AHSL(Unknown value)
        {
            return new AHSL(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(LAB value)
        {
            return new Unknown(value.L, value.A, value.B);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator LAB(Unknown value)
        {
            return new LAB(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(XYZ value)
        {
            return new Unknown(value.X, value.Y, value.Z);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator XYZ(Unknown value)
        {
            return new XYZ(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(CMYK value)
        {
            return new Unknown(value.Cyan, value.Magenta, value.Yellow);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator CMYK(Unknown value)
        {
            return new CMYK(value.X, value.Y, value.Z, 0);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(HSB value)
        {
            return new Unknown(value.Hue, value.Saturation, value.Brightness);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator HSB(Unknown value)
        {
            return new HSB(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(HSL value)
        {
            return new Unknown(value.Hue, value.Saturation, value.Lightness);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator HSL(Unknown value)
        {
            return new HSL(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(RGB value)
        {
            return new Unknown(value.Red, value.Green, value.Blue);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator RGB(Unknown value)
        {
            return new RGB(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(RYB value)
        {
            return new Unknown(value.Red, value.Yellow, value.Blue);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator RYB(Unknown value)
        {
            return new RYB(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(sRGB value)
        {
            return new Unknown(value.Red, value.Green, value.Blue);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator sRGB(Unknown value)
        {
            return new sRGB(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(YCbCr value)
        {
            return new Unknown(value.Y, value.Cb, value.Cr);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator YCbCr(Unknown value)
        {
            return new YCbCr(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(YCgCo value)
        {
            return new Unknown(value.Y, value.Cg, value.Co);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator YCgCo(Unknown value)
        {
            return new YCgCo(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(YDbDr value)
        {
            return new Unknown(value.Y, value.Db, value.Dr);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator YDbDr(Unknown value)
        {
            return new YDbDr(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(YIQ value)
        {
            return new Unknown(value.Y, value.I, value.Q);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator YIQ(Unknown value)
        {
            return new YIQ(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(YPbPr value)
        {
            return new Unknown(value.Y, value.Pb, value.Pr);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator YPbPr(Unknown value)
        {
            return new YPbPr(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(YUV value)
        {
            return new Unknown(value.Y, value.U, value.V);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator YUV(Unknown value)
        {
            return new YUV(value.X, value.Y, value.Z);
        }

        /// <summary>
        /// Defines an explicit conversion Space в Unknown.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Unknown(Color value)
        {
            return new Unknown(value.R, value.G, value.B);
        }
        /// <summary>
        /// Defines an explicit conversion Unknown в Space.
        /// </summary>
        /// <param name="value">Structure</param>
        /// <returns>Structure</returns>
        public static implicit operator Color(Unknown value)
        {
            return Color.FromArgb((int)value.X, (int)value.Y, (int)value.Z);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new Unknown(this.X, this.Y, this.Z);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public Unknown Clone()
        {
            return new Unknown(this.X, this.Y, this.Z);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Returns the color model RGB.
        /// </summary>
        public RGB ToRGB
        {
            get
            {
                return this;
            }
        }
        #endregion
    }
}
