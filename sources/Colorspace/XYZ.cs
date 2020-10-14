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
        public static readonly XYZ White = new XYZ(0.9505, 1.0, 1.0890);
        #endregion

        #region Private data
        private double x;
        private double y;
        private double z;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure CIE XYZ.
        /// </summary>
        /// <param name="x">Component X [0, 1]</param>
        /// <param name="y">Component Y [0, 1]</param>
        /// <param name="z">Component Z [0, 1]</param>
        public XYZ(double x, double y, double z)
        {
            this.x = (x > 1.0) ? 1.0 : ((x < 0) ? 0 : x);
            this.y = (y > 1.0) ? 1.0 : ((y < 0) ? 0 : y);
            this.z = (z > 1.0) ? 1.0 : ((z < 0) ? 0 : z);
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public double X
        {
            get
            {
                return this.x;
            }
            set
            {
                this.x = (value > 1.0) ? 1.0 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public double Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = (value > 1.0) ? 1.0 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public double Z
        {
            get
            {
                return this.z;
            }
            set
            {
                this.z = (value > 1.0) ? 1.0 : ((value < 0) ? 0 : value);
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
            return x.ToString() + "\n" + y.ToString() + "\n" + z.ToString();
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
        public static LAB ToLAB(double x, double y, double z)
        {
            LAB lab = new LAB();

            lab.L = 116.0 * Fxyz(y / XYZ.White.Y) - 16;
            lab.A = 500.0 * (Fxyz(x / XYZ.White.X) - Fxyz(y / XYZ.White.Y));
            lab.B = 200.0 * (Fxyz(y / XYZ.White.Y) - Fxyz(z / XYZ.White.Z));

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
            double rLinear = (double)red / 255.0;
            double gLinear = (double)green / 255.0;
            double bLinear = (double)blue / 255.0;

            // convert to a sRGB form
            double r = (rLinear > 0.04045) ? Math.Pow((rLinear + 0.055) / (1 + 0.055), 2.2) : (rLinear / 12.92);
            double g = (gLinear > 0.04045) ? Math.Pow((gLinear + 0.055) / (1 + 0.055), 2.2) : (gLinear / 12.92);
            double b = (bLinear > 0.04045) ? Math.Pow((bLinear + 0.055) / (1 + 0.055), 2.2) : (bLinear / 12.92);

            // converts
            return new XYZ(
                (r * 0.4124 + g * 0.3576 + b * 0.1805),
                (r * 0.2126 + g * 0.7152 + b * 0.0722),
                (r * 0.0193 + g * 0.1192 + b * 0.9505)
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
        /// 
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        private static double Fxyz(double t)
        {
            return ((t > 0.008856) ? Math.Pow(t, (1.0 / 3.0)) : (7.787 * t + 16.0 / 116.0));
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
                double[] linear = new double[3];
                linear[0] = x * 3.2410 - y * 1.5374 - z * 0.4986; // red
                linear[1] = -x * 0.9692 + y * 1.8760 - z * 0.0416; // green
                linear[2] = x * 0.0556 - y * 0.2040 + z * 1.0570; // blue
                double pow = 1.0 / 2.4;

                for (int i = 0; i < 3; i++)
                {
                    linear[i] = (linear[i] <= 0.0031308) ? 12.92 * linear[i] : (1 + 0.055) * Math.Pow(linear[i], pow) - 0.055;
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
