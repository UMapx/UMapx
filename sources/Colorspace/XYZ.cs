using System;

namespace UMapx.Colorspace
{
    /// <summary>
    /// Defines a color model CIE XYZ.
    /// </summary>
    [Serializable]
    public struct XYZ : IColorSpace, ICloneable
    {
        #region Readonly
        /// <summary>
        /// Returns white color.
        /// </summary>
        public static readonly XYZ White = new XYZ(0.9505f, 1.0f, 1.0890f);
        #endregion

        #region Private data
        private float x;
        private float y;
        private float z;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure CIE XYZ.
        /// </summary>
        /// <param name="x">Component X [0, 1]</param>
        /// <param name="y">Component Y [0, 1]</param>
        /// <param name="z">Component Z [0, 1]</param>
        public XYZ(float x, float y, float z)
        {
            this.x = (x > 1.0) ? 1.0f : ((x < 0) ? 0 : x);
            this.y = (y > 1.0) ? 1.0f : ((y < 0) ? 0 : y);
            this.z = (z > 1.0) ? 1.0f : ((z < 0) ? 0 : z);
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public float X
        {
            get
            {
                return this.x;
            }
            set
            {
                this.x = (value > 1.0) ? 1.0f : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public float Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = (value > 1.0) ? 1.0f : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public float Z
        {
            get
            {
                return this.z;
            }
            set
            {
                this.z = (value > 1.0) ? 1.0f : ((value < 0) ? 0 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">CIE XYZ structure</param>
        /// <param name="item2">CIE XYZ structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(XYZ item1, XYZ item2)
        {
            return (
                item1.x == item2.x
                && item1.y == item2.y
                && item1.z == item2.z
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">CIE XYZ structure</param>
        /// <param name="item2">CIE XYZ structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(XYZ item1, XYZ item2)
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

            return (this == (XYZ)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return x.GetHashCode() ^ y.GetHashCode() ^ z.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return $"{x}{Environment.NewLine}{y}{Environment.NewLine}{z}";
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new XYZ(this.X, this.Y, this.Z);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public XYZ Clone()
        {
            return new XYZ(this.X, this.Y, this.Z);
        }
        #endregion

        #region CIE XYZ convert
        /// <summary>
        /// Converts a color model CIE XYZ in model CIE Lab.
        /// </summary>
        /// <param name="x">Component X</param>
        /// <param name="y">Component Y</param>
        /// <param name="z">Component Z</param>
        /// <returns>CIE Lab structure</returns>
        public static LAB ToLAB(float x, float y, float z)
        {
            LAB lab = new LAB();

            lab.L = 116.0f * Fxyz(y / XYZ.White.Y) - 16;
            lab.A = 500.0f * (Fxyz(x / XYZ.White.X) - Fxyz(y / XYZ.White.Y));
            lab.B = 200.0f * (Fxyz(y / XYZ.White.Y) - Fxyz(z / XYZ.White.Z));

            return lab;
        }
        /// <summary>
        /// Converts a color model CIE XYZ in model CIE Lab.
        /// </summary>
        /// <param name="xyz">CIE XYZ structure</param>
        /// <returns>CIE Lab structure</returns>
        public static LAB ToLAB(XYZ xyz)
        {
            return ToLAB(xyz.X, xyz.Y, xyz.Z);
        }
        /// <summary>
        /// Converts a color model RGB in model CIE XYZ.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>CIE XYZ structure</returns>
        public static XYZ FromRGB(int red, int green, int blue)
        {
            // normalize red, green, blue values
            float rLinear = (float)red / 255.0f;
            float gLinear = (float)green / 255.0f;
            float bLinear = (float)blue / 255.0f;

            // convert to a sRGB form
            float r = (rLinear > 0.04045) ? (float)Math.Pow((rLinear + 0.055) / (1 + 0.055), 2.2) : (float)(rLinear / 12.92);
            float g = (gLinear > 0.04045) ? (float)Math.Pow((gLinear + 0.055) / (1 + 0.055), 2.2) : (float)(gLinear / 12.92);
            float b = (bLinear > 0.04045) ? (float)Math.Pow((bLinear + 0.055) / (1 + 0.055), 2.2) : (float)(bLinear / 12.92);

            // converts
            return new XYZ(
                (float)(r * 0.4124 + g * 0.3576 + b * 0.1805),
                (float)(r * 0.2126 + g * 0.7152 + b * 0.0722),
                (float)(r * 0.0193 + g * 0.1192 + b * 0.9505)
                );
        }
        /// <summary>
        /// Converts a color model RGB in model CIE XYZ.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>CIE XYZ structure</returns>
        public static XYZ FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Computes the nonlinear helper function used in CIE XYZ to LAB conversion.
        /// </summary>
        /// <param name="t">Input value</param>
        /// <returns>Transformed value</returns>
        private static float Fxyz(float t)
        {
            return ((t > 0.008856) ? (float)Math.Pow(t, (1.0 / 3.0)) : (float)(7.787 * t + 16.0 / 116.0));
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model CIE XYZ in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                float[] linear = new float[3];
                linear[0] = x * 3.2410f - y * 1.5374f - z * 0.4986f; // red
                linear[1] = -x * 0.9692f + y * 1.8760f - z * 0.0416f; // green
                linear[2] = x * 0.0556f - y * 0.2040f + z * 1.0570f; // blue
                float pow = 1.0f / 2.4f;

                for (int i = 0; i < 3; i++)
                {
                    linear[i] = (linear[i] <= 0.0031308) ? 12.92f * linear[i] : (1 + 0.055f) * (float)Math.Pow(linear[i], pow) - 0.055f;
                }

                return new RGB(
                    (int)(linear[0] * 255.0),
                    (int)(linear[1] * 255.0),
                    (int)(linear[2] * 255.0)
                    );
            }
        }
        #endregion
    }
}
