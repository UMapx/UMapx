using System;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model CIE Lab.
    /// </summary>
    [Serializable]
    public struct LAB : IColorSpace, ICloneable
    {
        #region Private data
        private float l;
        private float a;
        private float b;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure CIE Lab.
        /// </summary>
        /// <param name="l">Component L [0, 100]</param>
        /// <param name="a">Component a [-127, 127]</param>
        /// <param name="b">Component b [-127, 127]</param>
        public LAB(float l, float a, float b)
        {
            this.l = (l > 100.0) ? 100.0f : ((l < 0) ? 0 : l);
            this.a = (a > 127.0) ? 127.0f : ((a < -127) ? -127 : a);
            this.b = (b > 127.0) ? 127.0f : ((b < -127) ? -127 : b);
        }
        /// <summary>
        /// Defines a component of the model [0, 100].
        /// </summary>
        public float L
        {
            get
            {
                return this.l;
            }
            set
            {
                this.l = (value > 100.0) ? 100.0f : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [-127, 127].
        /// </summary>
        public float A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = (value > 127.0) ? 127.0f : ((value < -127) ? -127 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [-127, 127].
        /// </summary>
        public float B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = (value > 127.0) ? 127.0f : ((value < -127) ? -127 : value);
            }
        }
        #endregion

        #region Operators
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">CIE Lab structure</param>
        /// <param name="item2">CIE Lab structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(LAB item1, LAB item2)
        {
            return (
                item1.L == item2.L
                && item1.A == item2.A
                && item1.B == item2.B
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">CIE Lab structure</param>
        /// <param name="item2">CIE Lab structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(LAB item1, LAB item2)
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

            return (this == (LAB)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return l.GetHashCode() ^ a.GetHashCode() ^ b.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return $"{l}{Environment.NewLine}{a}{Environment.NewLine}{b}";
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new LAB(this.l, this.a, this.b);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public LAB Clone()
        {
            return new LAB(this.l, this.a, this.b);
        }
        #endregion

        #region CIE Lab convert
        /// <summary>
        /// Converts from CIE Lab to CIE XYZ.
        /// </summary>
        /// <param name="l">Component L</param>
        /// <param name="a">Component a</param>
        /// <param name="b">Component b</param>
        /// <returns>CIE XYZ structure</returns>
        public static XYZ ToXYZ(float l, float a, float b)
        {
            float theta = 6.0f / 29.0f;
            float fy = (l + 16) / 116.0f;
            float fx = fy + (a / 500.0f);
            float fz = fy - (b / 200.0f);
            float theta2 = theta * theta;
            float k = 16.0f / 116.0f;
            XYZ D65 = XYZ.White;

            return new XYZ(
                (fx > theta) ? D65.X * (fx * fx * fx) : (fx - k) * 3 * theta2 * D65.X,
                (fy > theta) ? D65.Y * (fy * fy * fy) : (fy - k) * 3 * theta2 * D65.Y,
                (fz > theta) ? D65.Z * (fz * fz * fz) : (fz - k) * 3 * theta2 * D65.Z
                );
        }
        /// <summary>
        /// Converts from CIE Lab to CIE XYZ.
        /// </summary>
        /// <param name="lab">CIE Lab structure</param>
        /// <returns>CIE XYZ structure</returns>
        public static XYZ ToXYZ(LAB lab)
        {
            return LAB.ToXYZ(lab.L, lab.A, lab.B);
        }
        /// <summary>
        /// Converts from RGB to CIE Lab.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>CIE Lab structure</returns>
        public static LAB FromRGB(int red, int green, int blue)
        {
            return XYZ.ToLAB(XYZ.FromRGB(red, green, blue));
        }
        /// <summary>
        /// Converts from RGB to CIE Lab.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>CIE Lab structure</returns>
        public static LAB FromRGB(RGB rgb)
        {
            return XYZ.ToLAB(XYZ.FromRGB(rgb.Red, rgb.Green, rgb.Blue));
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts from CIE Lab to RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                return LAB.ToXYZ(l, a, b).ToRGB;
            }
        }
        #endregion
    }
}
