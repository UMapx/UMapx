// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Drawing;
using UMapx.Core;

namespace UMapx.Colorspace
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                             UMAPX.COLORSPACE
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Colorspaces
    /// <summary>
    /// Defines a color model CIE Lab.
    /// </summary>
    [Serializable]
    public struct LAB : IColorSpace, ICloneable
    {
        #region Private data
        private double l;
        private double a;
        private double b;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure CIE Lab.
        /// </summary>
        /// <param name="l">Component L [0, 100]</param>
        /// <param name="a">Component a [-127, 127]</param>
        /// <param name="b">Component b [-127, 127]</param>
        public LAB(double l, double a, double b)
        {
            this.l = (l > 100.0) ? 100.0 : ((l < 0) ? 0 : l);
            this.a = (a > 127.0) ? 127.0 : ((a < -127) ? -127 : a);
            this.b = (b > 127.0) ? 127.0 : ((b < -127) ? -127 : b);
        }
        /// <summary>
        /// Defines a component of the model [0, 100].
        /// </summary>
        public double L
        {
            get
            {
                return this.l;
            }
            set
            {
                this.l = (value > 100.0) ? 100.0 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [-127, 127].
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = (value > 127.0) ? 127.0 : ((value < -127) ? -127 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [-127, 127].
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = (value > 127.0) ? 127.0 : ((value < -127) ? -127 : value);
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
            return L.GetHashCode() ^ a.GetHashCode() ^ b.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return l.ToString() + "\n" + a.ToString() + "\n" + b.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new LAB(this.L, this.A, this.B);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public LAB Clone()
        {
            return new LAB(this.L, this.A, this.B);
        }
        #endregion

        #region CIE Lab convert
        /// <summary>
        /// Converts a color model CIE Lab in model CIE XYZ.
        /// </summary>
        /// <param name="l">Component L</param>
        /// <param name="a">Component a</param>
        /// <param name="b">Component b</param>
        /// <returns>CIE XYZ structure</returns>
        public static XYZ ToXYZ(double l, double a, double b)
        {
            double theta = 6.0 / 29.0;
            double fy = (l + 16) / 116.0;
            double fx = fy + (a / 500.0);
            double fz = fy - (b / 200.0);
            double theta2 = theta * theta;
            double k = 16.0 / 116.0;
            XYZ D65 = XYZ.White;

            return new XYZ(
                (fx > theta) ? D65.X * (fx * fx * fx) : (fx - k) * 3 * theta2 * D65.X,
                (fy > theta) ? D65.Y * (fy * fy * fy) : (fy - k) * 3 * theta2 * D65.Y,
                (fz > theta) ? D65.Z * (fz * fz * fz) : (fz - k) * 3 * theta2 * D65.Z
                );
        }
        /// <summary>
        /// Converts a color model CIE Lab in model CIE XYZ.
        /// </summary>
        /// <param name="lab">CIE Lab structure</param>
        /// <returns>CIE XYZ structure</returns>
        public static XYZ ToXYZ(LAB lab)
        {
            return LAB.ToXYZ(lab.L, lab.A, lab.B);
        }
        /// <summary>
        /// Converts a color model RGB in model CIE Lab.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>CIE Lab structure</returns>
        public static LAB ToLAB(int red, int green, int blue)
        {
            return XYZ.ToLAB(XYZ.FromRGB(red, green, blue));
        }
        /// <summary>
        /// Converts a color model RGB in model CIE Lab.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>CIE Lab structure</returns>
        public static LAB ToLAB(RGB rgb)
        {
            return XYZ.ToLAB(XYZ.FromRGB(rgb.Red, rgb.Green, rgb.Blue));
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model CIE Lab in model RGB.
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

            lab.L = 116.0 *  Fxyz(y / XYZ.White.Y) - 16;
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
    /// <summary>
    /// Defines a color model СMYK.
    /// </summary>
    [Serializable]
    public struct CMYK : IColorSpace, ICloneable
    {
        #region Private data
        private double c;
        private double m;
        private double y;
        private double k;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure CMYK.
        /// </summary>
        /// <param name="c">Cyan [0, 1]</param>
        /// <param name="m">Magenta [0, 1]</param>
        /// <param name="y">Yellow [0, 1]</param>
        /// <param name="k">Keycolor [0, 1]</param>
        public CMYK(double c, double m, double y, double k)
        {
            this.c = (c > 1) ? 1 : ((c < 0) ? 0 : c);
            this.m = (m > 1) ? 1 : ((m < 0) ? 0 : m);
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.k = (k > 1) ? 1 : ((k < 0) ? 0 : k);
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public double Cyan
        {
            get
            {
                return c;
            }
            set
            {
                c = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public double Magenta
        {
            get
            {
                return m;
            }
            set
            {
                m = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public double Yellow
        {
            get
            {
                return y;
            }
            set
            {
                y = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the model [0, 1].
        /// </summary>
        public double Keycolor
        {
            get
            {
                return k;
            }
            set
            {
                k = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">CMYK structure</param>
        /// <param name="item2">CMYK structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(CMYK item1, CMYK item2)
        {
            return (
                item1.Cyan == item2.Cyan
                && item1.Magenta == item2.Magenta
                && item1.Yellow == item2.Yellow
                && item1.Keycolor == item2.Keycolor
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">CMYK structure</param>
        /// <param name="item2">CMYK structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(CMYK item1, CMYK item2)
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

            return (this == (CMYK)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Cyan.GetHashCode() ^ Magenta.GetHashCode() ^ Yellow.GetHashCode() ^ Keycolor.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Cyan.ToString() + "\n" + Magenta.ToString() + "\n" + Yellow.ToString() + "\n" + Keycolor.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new CMYK(this.Cyan, this.Magenta, this.Yellow, this.Keycolor);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public CMYK Clone()
        {
            return new CMYK(this.Cyan, this.Magenta, this.Yellow, this.Keycolor);
        }
        #endregion

        #region CMYK convert
        /// <summary>
        /// Converts a color model RGB in model CMYK.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>CMYK structure</returns>
        public static CMYK FromRGB(int red, int green, int blue)
        {
            double c = (255.0f - red) / 255.0f;
            double m = (255.0f - green) / 255.0f;
            double y = (255.0f - blue) / 255.0f;

            double min = Maths.Min(c, m, y);

            if (min == 1.0)
            {
                return new CMYK(0.0f, 0.0f, 0.0f, 1.0f);
            }
            else
            {
                double k = 1.0f - min;
                return new CMYK((c - min) / k, (m - min) / k, (y - min) / k, min);
            }
        }
        /// <summary>
        /// Converts a color model RGB in model HSB.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>CMYK structure</returns>
        public static CMYK FromRGB(RGB rgb)
        {
            return CMYK.FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model CMYK in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                double k1 = (1.0f - k) * 255.0;

                int r = (int)((1.0f - c) * k1);
                int g = (int)((1.0f - m) * k1);
                int b = (int)((1.0f - y) * k1);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model HSB.
    /// </summary>
    [Serializable]
    public struct HSB : IColorSpace, ICloneable
    {
        #region Private data
        private double h;
        private double s;
        private double b;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure HSB.
        /// </summary>
        /// <param name="h">Hue [0, 359]</param>
        /// <param name="s">Saturation [0, 1]</param>
        /// <param name="b">Brightness [0, 1]</param>
        public HSB(double h, double s, double b)
        {
            this.h = (h > 359) ? 359 : ((h < 0) ? 0 : h);
            this.s = (s > 1) ? 1 : ((s < 0) ? 0 : s);
            this.b = (b > 1) ? 1 : ((b < 0) ? 0 : b);
        }
        /// <summary>
        /// Defines a component of the color model [0, 359].
        /// </summary>
        public double Hue
        {
            get
            {
                return h;
            }
            set
            {
                h = (value > 359) ? 359 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Saturation
        {
            get
            {
                return s;
            }
            set
            {
                s = (value > 1) ? 1 : ((value < 0.00001f) ? 0.00001f : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Brightness
        {
            get
            {
                return b;
            }
            set
            {
                b = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">HSB structure</param>
        /// <param name="item2">HSB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(HSB item1, HSB item2)
        {
            return (
                item1.Hue == item2.Hue
                && item1.Saturation == item2.Saturation
                && item1.Brightness == item2.Brightness
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">HSB structure</param>
        /// <param name="item2">HSB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(HSB item1, HSB item2)
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

            return (this == (HSB)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Hue.GetHashCode() ^ Saturation.GetHashCode() ^ Brightness.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Hue.ToString() + "\n" + Saturation.ToString() + "\n" + Brightness.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new HSB(this.Hue, this.Saturation, this.Brightness);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public HSB Clone()
        {
            return new HSB(this.Hue, this.Saturation, this.Brightness);
        }
        #endregion

        #region HSB convert
        /// <summary>
        /// Converts a color model RGB in model HSB.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>HSB structure</returns>
        public static HSB FromRGB(int red, int green, int blue)
        {
            double r = red / 255.0f;
            double g = green / 255.0f;
            double b = blue / 255.0f;

            double max = Maths.Max(r, g, b);
            double min = Maths.Min(r, g, b);

            double h = 0;
            double l = max - min;

            if (max == r && g >= b)
            {
                if (l == 0.0) h = 0;
                else h = (int)(60 * (g - b) / l);
            }
            else if (max == r && g < b)
            {
                h = (int)(60 * (g - b) / l + 360);
            }
            else if (max == g)
            {
                h = (int)(60 * (b - r) / l + 120);
            }
            else if (max == b)
            {
                h = (int)(60 * (r - g) / l + 240);
            }

            double s = (max == 0.0f) ? 0.0f : (1 - (min / max));
            return new HSB(h, s, max);
        }
        /// <summary>
        /// Converts a color model RGB in model HSB.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>HSB structure</returns>
        public static HSB FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model HSB in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                double red = 0;
                double green = 0;
                double blue = 0;

                if (b <= 0)
                {
                    red = green = blue = 0;
                }

                if (s <= 0)
                {
                    red = green = blue = b;
                }
                else
                {
                    double Hi = h / 60.0f;
                    int HiMod6 = (int)(Math.Floor(Hi));
                    double Sector = Hi - HiMod6;

                    double Vmin = b * (1 - s);
                    double Vdec = b * (1 - (s * Sector));
                    double Vinc = b * (1 - (s * (1 - Sector)));

                    switch (HiMod6)
                    {
                        case 0:
                            red = b;
                            green = Vinc;
                            blue = Vmin;
                            break;
                        case 1:
                            red = Vdec;
                            green = b;
                            blue = Vmin;
                            break;
                        case 2:
                            red = Vmin;
                            green = b;
                            blue = Vinc;
                            break;
                        case 3:
                            red = Vmin;
                            green = Vdec;
                            blue = b;
                            break;
                        case 4:
                            red = Vinc;
                            green = Vmin;
                            blue = b;
                            break;
                        case 5:
                            red = b;
                            green = Vmin;
                            blue = Vdec;
                            break;
                    }
                }
                return new RGB(red * 255, green * 255, blue * 255);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model HSL.
    /// </summary>
    [Serializable]
    public struct HSL : IColorSpace, ICloneable
    {
        #region Private data
        private double h;
        private double s;
        private double l;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure HSL.
        /// </summary>
        /// <param name="h">Hue [0, 360]</param>
        /// <param name="s">Saturation [0, 1]</param>
        /// <param name="l">Lightness [0, 1]</param>
        public HSL(double h, double s, double l)
        {
            this.h = (h > 360) ? 360 : ((h < 0) ? 0 : h);
            this.s = (s > 1) ? 1 : ((s < 0) ? 0 : s);
            this.l = (l > 1) ? 1 : ((l < 0) ? 0 : l);
        }
        /// <summary>
        /// Defines a component of the color model [0, 360].
        /// </summary>
        public double Hue
        {
            get
            {
                return h;
            }
            set
            {
                h = (value > 360) ? 360 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Saturation
        {
            get
            {
                return s;
            }
            set
            {
                s = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Lightness
        {
            get
            {
                return l;
            }
            set
            {
                l = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">HSL structure</param>
        /// <param name="item2">HSL structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(HSL item1, HSL item2)
        {
            return (
                item1.Hue == item2.Hue
                && item1.Saturation == item2.Saturation
                && item1.Lightness == item2.Lightness
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">HSL structure</param>
        /// <param name="item2">HSL structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(HSL item1, HSL item2)
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

            return (this == (HSL)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Hue.GetHashCode() ^ Saturation.GetHashCode() ^ Lightness.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Hue.ToString() + "\n" + Saturation.ToString() + "\n" + Lightness.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new HSL(this.Hue, this.Saturation, this.Lightness);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public HSL Clone()
        {
            return new HSL(this.Hue, this.Saturation, this.Lightness);
        }
        #endregion

        #region HSL convert
        /// <summary>
        /// Converts a color model RGB in model HSL.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>HSL structure</returns>
        public static HSL FromRGB(int red, int green, int blue)
        {
            double s = 0.0f, l = 0.0f;
            int h = 0;

            double r = red / 255.0f;
            double g = green / 255.0f;
            double b = blue / 255.0f;

            double max = Maths.Max(r, g, b);
            double min = Maths.Min(r, g, b);

            if (max == min)
            {
                h = 0;
            }
            else if (max == r && g >= b)
            {
                h = (int)(60 * (g - b) / (max - min));
            }
            else if (max == r && g < b)
            {
                h = (int)(60 * (g - b) / (max - min) + 360);
            }
            else if (max == g)
            {
                h = (int)(60 * (b - r) / (max - min) + 120);
            }
            else if (max == b)
            {
                h = (int)(60 * (r - g) / (max - min) + 240);
            }

            l = (max + min) / 2.0f;

            if (l == 0.0 || max == min)
            {
                s = 0.0f;
            }
            else if (0.0 < l && l <= 0.5)
            {
                s = (max - min) / (max + min);
            }
            else if (l > 0.5)
            {
                s = (max - min) / (2.0f - (max + min));
            }

            return new HSL(h, s, l);
        }
        /// <summary>
        /// Converts a color model RGB in model HSL.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>HSL structure</returns>
        public static HSL FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model HSL in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                if (s == 0)
                {
                    double rgb = l * 255;
                    return new RGB(rgb, rgb, rgb);
                }
                else
                {
                    double q = (l < 0.5f) ? (l * (1.0f + s)) : (l + s - (l * s));
                    double p = (2.0f * l) - q;

                    double Hk = h / 360.0f;
                    double[] T = new double[3];

                    T[0] = Hk + 0.3333f;	// Tr
                    T[1] = Hk;				// Tb
                    T[2] = Hk - 0.3333f;	// Tg

                    for (int i = 0; i < 3; i++)
                    {
                        if (T[i] < 0) T[i] += 1.0f;
                        if (T[i] > 1) T[i] -= 1.0f;

                        if ((T[i] * 6) < 1)
                        {
                            T[i] = p + ((q - p) * 6.0f * T[i]);
                        }
                        else if ((T[i] * 2.0) < 1)
                        {
                            T[i] = q;
                        }
                        else if ((T[i] * 3.0) < 2)
                        {
                            T[i] = p + (q - p) * ((2.0f / 3.0f) - T[i]) * 6.0f;
                        }
                        else T[i] = p;
                    }

                    double red = T[0] * 255;
                    double green = T[1] * 255;
                    double blue = T[2] * 255;

                    return new RGB(red, green, blue);
                }
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model AHSL.
    /// </summary>
    [Serializable]
    public struct AHSL : IColorSpace, ICloneable
    {
        #region Private data
        private double h;
        private double s;
        private double l;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure HSL.
        /// </summary>
        /// <param name="h">Hue [0, 360]</param>
        /// <param name="s">Saturation [0, 255]</param>
        /// <param name="l">Lightness [-100, 100]</param>
        public AHSL(double h, double s, double l)
        {
            this.h = (h > 359) ? 359 : ((h < 0) ? 0 : h);
            this.s = (s > 255) ? 255 : ((s < 0) ? 0 : s);
            this.l = (l > 100) ? 100 : ((l < -100) ? -100 : l);
        }
        /// <summary>
        /// Defines a component of the color model [0, 359].
        /// </summary>
        public double Hue
        {
            get
            {
                return h;
            }
            set
            {
                h = (value > 359) ? 359 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 255].
        /// </summary>
        public double Saturation
        {
            get
            {
                return s;
            }
            set
            {
                s = (value > 255) ? 255 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-100, 100].
        /// </summary>
        public double Lightness
        {
            get
            {
                return l;
            }
            set
            {
                l = (value > 100) ? 100 : ((value < -100) ? -100 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">HSL structure</param>
        /// <param name="item2">HSL structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(AHSL item1, AHSL item2)
        {
            return (
                item1.Hue == item2.Hue
                && item1.Saturation == item2.Saturation
                && item1.Lightness == item2.Lightness
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">HSL structure</param>
        /// <param name="item2">HSL structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(AHSL item1, AHSL item2)
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

            return (this == (AHSL)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Hue.GetHashCode() ^ Saturation.GetHashCode() ^ Lightness.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Hue.ToString() + "\n" + Saturation.ToString() + "\n" + Lightness.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new AHSL(this.Hue, this.Saturation, this.Lightness);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public AHSL Clone()
        {
            return new AHSL(this.Hue, this.Saturation, this.Lightness);
        }
        #endregion

        #region AHSL covnert
        /// <summary>
        /// Converts a color model RGB in model AHSL.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>HSL structure</returns>
        public static AHSL FromRGB(int red, int green, int blue)
        {
            double s = 0, l = 0, h = 0;

            double max = Maths.Max(red, green, blue);
            double min = Maths.Min(red, green, blue);

            if (max == min)
            {
                h = 0;
            }
            else if (max == red && green >= blue)
            {
                h = (int)(60 * (green - blue) / (max - min));
            }
            else if (max == red && green < blue)
            {
                h = (int)(60 * (green - blue) / (max - min) + 360);
            }
            else if (max == green)
            {
                h = (int)(60 * (blue - red) / (max - min) + 120);
            }
            else if (max == blue)
            {
                h = (int)(60 * (red - green) / (max - min) + 240);
            }

            double gray = (red + green + blue) / 3.0;
            double r0 = 0, g0 = 0, b0 = 0;

            if ((h >= 0) && (h <= 60)) { r0 = 255.0; g0 = 4.25 * h; }
            else if ((h > 60) && (h <= 120)) { g0 = 255.0; r0 = 255 - 4.25 * (h - 60); }
            else if ((h > 120) && (h <= 180)) { g0 = 255.0; b0 = 4.25 * (h - 120); }
            else if ((h > 180) && (h <= 240)) { b0 = 255.0; g0 = 255 - 4.25 * (h - 180); }
            else if ((h > 240) && (h <= 300)) { b0 = 255.0; r0 = 4.25 * (h - 240); }
            else if ((h > 300) && (h <= 360)) { r0 = 255.0; b0 = 255 - 4.25 * (h - 300); }

            double gray0 = (r0 + g0 + b0) / 3.0;

            if (gray == gray0) { l = 0; }
            else if (gray > gray0) { l = 100 * (gray - gray0) / (255.0 - gray0); }
            else if (gray < gray0) { l = 100 * (gray - gray0) / gray0; }

            if (l > 0) { r0 = r0 + l * (255 - r0) / 100.0; }
            else if (l < 0) { r0 = r0 + l * r0 / 100.0; }

            if (red == gray) { s = 0; }
            else { s = 255 * Math.Abs(red - gray) / (Math.Abs(r0 - gray)); }

            return new AHSL(h, s, l);
        }
        /// <summary>
        /// Converts a color model RGB in model AHSL.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>HSL structure</returns>
        public static AHSL FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model AHSL in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                double r = 0, g = 0, b = 0;

                if ((h >= 0) && (h <= 60)) { r = 255.0; g = 4.25 * h; }
                else if ((h > 60) && (h <= 120)) { g = 255.0; r = 255 - 4.25 * (h - 60); }
                else if ((h > 120) && (h <= 180)) { g = 255.0; b = 4.25 * (h - 120); }
                else if ((h > 180) && (h <= 240)) { b = 255.0; g = 255 - 4.25 * (h - 180); }
                else if ((h > 240) && (h <= 300)) { b = 255.0; r = 4.25 * (h - 240); }
                else if ((h > 300) && (h <= 360)) { r = 255.0; b = 255 - 4.25 * (h - 300); }

                if (l > 0)
                {
                    r = r + l * (255 - r) / 100;
                    g = g + l * (255 - g) / 100;
                    b = b + l * (255 - b) / 100;
                }
                else
                {
                    r = r + l * r / 100.0;
                    g = g + l * g / 100.0;
                    b = b + l * b / 100.0;
                }

                double gray = (r + g + b) / 3.0;
                r = gray + s * (r - gray) / 255.0;
                g = gray + s * (g - gray) / 255.0;
                b = gray + s * (b - gray) / 255.0;

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model RGB.
    /// </summary>
    [Serializable]
    public struct RGB : IColorSpace, ICloneable
    {
        #region Private data
        private byte r;
        private byte g;
        private byte b;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure RGB.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        public RGB(int red, int green, int blue)
        {
            this.r = (byte)((red > 255) ? 255 : ((red < 0) ? 0 : red));
            this.g = (byte)((green > 255) ? 255 : ((green < 0) ? 0 : green));
            this.b = (byte)((blue > 255) ? 255 : ((blue < 0) ? 0 : blue));
        }
        /// <summary>
        /// Creates an instance of the structure RGB.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        public RGB(double red, double green, double blue)
        {
            this.r = (byte)((red > 255) ? 255 : ((red < 0) ? 0 : red));
            this.g = (byte)((green > 255) ? 255 : ((green < 0) ? 0 : green));
            this.b = (byte)((blue > 255) ? 255 : ((blue < 0) ? 0 : blue));
        }
        /// <summary>
        /// Defines a component of the color model [0, 255].
        /// </summary>
        public byte Red
        {
            get
            {
                return r;
            }
            set
            {
                r = (byte)((value > 255) ? 255 : ((value < 0) ? 0 : value));
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 255].
        /// </summary>
        public byte Green
        {
            get
            {
                return g;
            }
            set
            {
                g = (byte)((value > 255) ? 255 : ((value < 0) ? 0 : value));
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 255].
        /// </summary>
        public byte Blue
        {
            get
            {
                return b;
            }
            set
            {
                b = (byte)((value > 255) ? 255 : ((value < 0) ? 0 : value));
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">RGB structure</param>
        /// <param name="item2">RGB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(RGB item1, RGB item2)
        {
            return (
                item1.Red == item2.Red
                && item1.Green == item2.Green
                && item1.Blue == item2.Blue
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">RGB structure</param>
        /// <param name="item2">RGB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(RGB item1, RGB item2)
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

            return (this == (RGB)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Red.GetHashCode() ^ Green.GetHashCode() ^ Blue.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Red.ToString() + "\n" + Green.ToString() + "\n" + Blue.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new RGB(this.Red, this.Green, this.Blue);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public RGB Clone()
        {
            return new RGB(this.Red, this.Green, this.Blue);
        }
        #endregion

        #region Conversion operators
        /// <summary>
        /// Defines an explicit conversion RGB в System.Drawing.Color.
        /// </summary>
        /// <param name="value">RGB structure</param>
        /// <returns>Color in terms of red, green and blue</returns>
        public static implicit operator Color(RGB value)
        {
            return Color.FromArgb(value.Red, value.Green, value.Blue);
        }
        /// <summary>
        /// Defines an explicit conversion RGB в System.Drawing.Color.
        /// </summary>
        /// <param name="value">Color in terms of red, green and blue</param>
        /// <returns>RGB structure</returns>
        public static implicit operator RGB(Color value)
        {
            return new RGB(value.R, value.G, value.B);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Gets the int equivalent for a hexadecimal value.
        /// </summary>
        private static int GetIntFromHex(string strHex)
        {
            switch (strHex)
            {
                case ("A"):
                    {
                        return 10;
                    }
                case ("B"):
                    {
                        return 11;
                    }
                case ("C"):
                    {
                        return 12;
                    }
                case ("D"):
                    {
                        return 13;
                    }
                case ("E"):
                    {
                        return 14;
                    }
                case ("F"):
                    {
                        return 15;
                    }
                default:
                    {
                        return int.Parse(strHex);
                    }
            }
        }
        #endregion

        #region HEX convert
        /// <summary>
        /// Converts a color model HEX in model RGB.
        /// </summary>
        /// <param name="hexColor">HEX</param>
        /// <returns>RGB structure</returns>
        public static RGB FromHEX(string hexColor)
        {
            string r, g, b;

            hexColor = hexColor.Trim();
            if (hexColor[0] == '#') hexColor = hexColor.Substring(1, hexColor.Length - 1);

            r = hexColor.Substring(0, 2);
            g = hexColor.Substring(2, 2);
            b = hexColor.Substring(4, 2);

            r = Convert.ToString(16 * GetIntFromHex(r.Substring(0, 1)) + GetIntFromHex(r.Substring(1, 1)));
            g = Convert.ToString(16 * GetIntFromHex(g.Substring(0, 1)) + GetIntFromHex(g.Substring(1, 1)));
            b = Convert.ToString(16 * GetIntFromHex(b.Substring(0, 1)) + GetIntFromHex(b.Substring(1, 1)));

            return new RGB(
                Convert.ToInt32(r),
                Convert.ToInt32(g),
                Convert.ToInt32(b)
                );
        }
        /// <summary>
        /// Converts a color model RGB in model HEX.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string ToHEX(int red, int green, int blue)
        {
            return String.Format("#{0:x2}{1:x2}{2:x2}", red, green, blue);
        }
        /// <summary>
        /// Converts a color model RGB in model HEX.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string ToHEX(RGB rgb)
        {
            return ToHEX(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region Average
        /// <summary>
        /// Calculates the average brightness value.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>Double precision floating point number</returns>
        public static int Average(int red, int green, int blue)
        {
            return (red + green + blue) / 3;
        }
        /// <summary>
        /// Calculates the average brightness value.
        /// </summary>
        /// <param name="red">Red</param>
        /// <param name="green">Green</param>
        /// <param name="blue">Blue</param>
        /// <returns>Double precision floating point number</returns>
        public static double Average(double red, double green, double blue)
        {
            return (red + green + blue) / 3.0;
        }
        /// <summary>
        /// Calculates the average brightness value.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>Double precision floating point number</returns>
        public static int Average(RGB rgb)
        {
            return Average(rgb.Red, rgb.Green, rgb.Blue);
        }
        /// <summary>
        /// Calculates the average brightness value.
        /// </summary>
        /// <param name="rgb">sRGB structure</param>
        /// <returns>Double precision floating point number</returns>
        public static double Average(sRGB rgb)
        {
            return Average(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region PAL/NTC
        /// <summary>
        /// Calculates the brightness value in the standard (PAL/NTC).
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>Double precision floating point number</returns>
        public static int PAL(int red, int green, int blue)
        {
            return (int)(0.299 * red + 0.587 * green + 0.114 * blue);
        }
        /// <summary>
        /// Calculates the brightness value in the standard (PAL/NTC).
        /// </summary>
        /// <param name="red">Red</param>
        /// <param name="green">Green</param>
        /// <param name="blue">Blue</param>
        /// <returns>Double precision floating point number</returns>
        public static double PAL(double red, double green, double blue)
        {
            return 0.299 * red + 0.587 * green + 0.114 * blue;
        }
        /// <summary>
        /// Calculates the brightness value in the standard (PAL/NTC).
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>Double precision floating point number</returns>
        public static double PAL(RGB rgb)
        {
            return PAL(rgb.Red, rgb.Green, rgb.Blue);
        }
        /// <summary>
        /// Calculates the brightness value in the standard (PAL/NTC).
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>Double precision floating point number</returns>
        public static double PAL(sRGB rgb)
        {
            return PAL(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region HDTV
        /// <summary>
        /// Calculates the brightness value in the standard HDTV.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>Double precision floating point number</returns>
        public static int HDTV(int red, int green, int blue)
        {
            return (int)(0.2126 * red + 0.7152 * green + 0.0722 * blue);
        }
        /// <summary>
        /// Calculates the brightness value in the standard HDTV.
        /// </summary>
        /// <param name="red">Red</param>
        /// <param name="green">Green</param>
        /// <param name="blue">Blue</param>
        /// <returns>Double precision floating point number</returns>
        public static double HDTV(double red, double green, double blue)
        {
            return 0.2126 * red + 0.7152 * green + 0.0722 * blue;
        }
        /// <summary>
        /// Calculates the brightness value in the standard HDTV.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>Double precision floating point number</returns>
        public static int HDTV(RGB rgb)
        {
            return HDTV(rgb.Red, rgb.Green, rgb.Blue);
        }
        /// <summary>
        /// Calculates the brightness value in the standard HDTV.
        /// </summary>
        /// <param name="rgb">sRGB structure</param>
        /// <returns>Double precision floating point number</returns>
        public static double HDTV(sRGB rgb)
        {
            return HDTV(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RYY
        /// <summary>
        /// Calculates the brightness value in the standard RYY.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>Double precision floating point number</returns>
        public static int RYY(int red, int green, int blue)
        {
            return (int)(0.5 * red + 0.419 * green + 0.081 * blue);
        }
        /// <summary>
        /// Calculates the brightness value in the standard RYY.
        /// </summary>
        /// <param name="red">Red</param>
        /// <param name="green">Green</param>
        /// <param name="blue">Blue</param>
        /// <returns>Double precision floating point number</returns>
        public static double RYY(double red, double green, double blue)
        {
            return 0.5 * red + 0.419 * green + 0.081 * blue;
        }
        /// <summary>
        /// Calculates the brightness value in the standard RYY.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>Double precision floating point number</returns>
        public static int RYY(RGB rgb)
        {
            return RYY(rgb.Red, rgb.Green, rgb.Blue);
        }
        /// <summary>
        /// Calculates the brightness value in the standard RYY.
        /// </summary>
        /// <param name="rgb">sRGB structure</param>
        /// <returns>Double precision floating point number</returns>
        public static double RYY(sRGB rgb)
        {
            return RYY(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region Temperature
        /// <summary>
        /// Converts temperature T (in kelvins) to color in terms of red, green, and blue channels.
        /// </summary>
        /// <param name="temperature">Temperature [1000K, 10000K]</param>
        /// <returns>RGB structure</returns>
        public static RGB Temp2RGB(double temperature)
        {
            // Approximation of Planckian locus in RGB model
            // Designed by Asiryan Valeriy, 2018

            double r = 0, g, b = 0;          // color channels,
            double x = temperature / 1000.0; // normalize temerature,
            double y = x * x, z = y * x;     // variables.

            // Approximate red channel:
            if (x < 5)
            {
                r = 255;
            }
            else if (x <= 7)
            {
                r = -1575 + 718.5 * x - 70.5 * y;
            }
            // ************************************************************************

            // Approximate green channel:
            g = -50.1667 + 99.9846 * x - 7.7844 * y - 0.1175 * z;
            // ************************************************************************

            // Aproximate blue channel:
            if (x > 6)
            {
                b = 255;
            }
            else if (x >= 3)
            {
                b = 2247 - 1706 * x + 409 * y - 30 * z;
            }
            // ************************************************************************

            // Result color:
            return new RGB(r, g, b);
        }
        #endregion

        #region Saturation
        /// <summary>
        /// Corrects color saturation.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <param name="s">Saturation</param>
        /// <returns>RGB structure</returns>
        public static RGB Saturation(int red, int green, int blue, double s)
        {
            double max = Maths.Max(red, green, blue);

            // Result color:
            return new RGB(
                Maths.Byte(blue + (blue - max) * s),
                Maths.Byte(green + (green - max) * s),
                Maths.Byte(red + (red - max) * s));
        }
        /// <summary>
        /// Corrects color saturation.
        /// </summary>
        /// <param name="red">Red </param>
        /// <param name="green">Green</param>
        /// <param name="blue">Blue</param>
        /// <param name="s">Saturation</param>
        /// <returns>RGB structure</returns>
        public static sRGB Saturation(double red, double green, double blue, double s)
        {
            double max = Maths.Max(red, green, blue);

            // Result color:
            return new sRGB(
                Maths.Byte(blue + (blue - max) * s),
                Maths.Byte(green + (green - max) * s),
                Maths.Byte(red + (red - max) * s));
        }
        /// <summary>
        /// Corrects color saturation.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <param name="s">Saturation</param>
        /// <returns>RGB structure</returns>
        public static RGB Saturation(RGB rgb, double s)
        {
            return Saturation(rgb.Red, rgb.Green, rgb.Blue, s);
        }
        /// <summary>
        /// Corrects color saturation.
        /// </summary>
        /// <param name="rgb">sRGB structure</param>
        /// <param name="s">Saturation</param>
        /// <returns>RGB structure</returns>
        public static sRGB Saturation(sRGB rgb, double s)
        {
            return Saturation(rgb.Red, rgb.Green, rgb.Blue, s);
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

        #region Scheme
        /// <summary>
        /// Generates a color scheme.
        /// </summary>
        /// <param name="hue">Hue [0, 360]</param>
        /// <param name="length">Length</param>
        /// <returns>Color scheme</returns>
        public static RGB[] SchemeFromHue(double hue, uint length)
        {
            RGB[] scheme = new RGB[length];
            hue %= 360.0;

            for (int i = 0; i < scheme.Length; i++)
            {
                scheme[i] = (new HSL(hue, 1.0, i / (scheme.Length - 1.0))).ToRGB;
            }

            return scheme;
        }
        #endregion

        #region Completed shemes
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] Cool
        {
            get
            {
                return new RGB[] { new RGB(0, 255, 255), new RGB(4, 251, 255), new RGB(8, 247, 255), new RGB(12, 243, 255), new RGB(16, 239, 255), new RGB(20, 235, 255), new RGB(24, 231, 255), new RGB(28, 227, 255), new RGB(32, 223, 255), new RGB(36, 219, 255), new RGB(40, 215, 255), new RGB(45, 210, 255), new RGB(49, 206, 255), new RGB(53, 202, 255), new RGB(57, 198, 255), new RGB(61, 194, 255), new RGB(65, 190, 255), new RGB(69, 186, 255), new RGB(73, 182, 255), new RGB(77, 178, 255), new RGB(81, 174, 255), new RGB(85, 170, 255), new RGB(89, 166, 255), new RGB(93, 162, 255), new RGB(97, 158, 255), new RGB(101, 154, 255), new RGB(105, 150, 255), new RGB(109, 146, 255), new RGB(113, 142, 255), new RGB(117, 138, 255), new RGB(121, 134, 255), new RGB(125, 130, 255), new RGB(130, 125, 255), new RGB(134, 121, 255), new RGB(138, 117, 255), new RGB(142, 113, 255), new RGB(146, 109, 255), new RGB(150, 105, 255), new RGB(154, 101, 255), new RGB(158, 97, 255), new RGB(162, 93, 255), new RGB(166, 89, 255), new RGB(170, 85, 255), new RGB(174, 81, 255), new RGB(178, 77, 255), new RGB(182, 73, 255), new RGB(186, 69, 255), new RGB(190, 65, 255), new RGB(194, 61, 255), new RGB(198, 57, 255), new RGB(202, 53, 255), new RGB(206, 49, 255), new RGB(210, 45, 255), new RGB(215, 40, 255), new RGB(219, 36, 255), new RGB(223, 32, 255), new RGB(227, 28, 255), new RGB(231, 24, 255), new RGB(235, 20, 255), new RGB(239, 16, 255), new RGB(243, 12, 255), new RGB(247, 8, 255), new RGB(251, 4, 255), new RGB(255, 0, 255) };
            }
        }
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] Hot
        {
            get
            {
                return new RGB[] { new RGB(11, 0, 0), new RGB(21, 0, 0), new RGB(32, 0, 0), new RGB(43, 0, 0), new RGB(53, 0, 0), new RGB(64, 0, 0), new RGB(74, 0, 0), new RGB(85, 0, 0), new RGB(96, 0, 0), new RGB(106, 0, 0), new RGB(117, 0, 0), new RGB(128, 0, 0), new RGB(138, 0, 0), new RGB(149, 0, 0), new RGB(159, 0, 0), new RGB(170, 0, 0), new RGB(181, 0, 0), new RGB(191, 0, 0), new RGB(202, 0, 0), new RGB(213, 0, 0), new RGB(223, 0, 0), new RGB(234, 0, 0), new RGB(244, 0, 0), new RGB(255, 0, 0), new RGB(255, 11, 0), new RGB(255, 21, 0), new RGB(255, 32, 0), new RGB(255, 43, 0), new RGB(255, 53, 0), new RGB(255, 64, 0), new RGB(255, 74, 0), new RGB(255, 85, 0), new RGB(255, 96, 0), new RGB(255, 106, 0), new RGB(255, 117, 0), new RGB(255, 128, 0), new RGB(255, 138, 0), new RGB(255, 149, 0), new RGB(255, 159, 0), new RGB(255, 170, 0), new RGB(255, 181, 0), new RGB(255, 191, 0), new RGB(255, 202, 0), new RGB(255, 213, 0), new RGB(255, 223, 0), new RGB(255, 234, 0), new RGB(255, 244, 0), new RGB(255, 255, 0), new RGB(255, 255, 16), new RGB(255, 255, 32), new RGB(255, 255, 48), new RGB(255, 255, 64), new RGB(255, 255, 80), new RGB(255, 255, 96), new RGB(255, 255, 112), new RGB(255, 255, 128), new RGB(255, 255, 143), new RGB(255, 255, 159), new RGB(255, 255, 175), new RGB(255, 255, 191), new RGB(255, 255, 207), new RGB(255, 255, 223), new RGB(255, 255, 239), new RGB(255, 255, 255) };
            }
        }
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] Copper
        {
            get
            {
                return new RGB[] { new RGB(0, 0, 0), new RGB(5, 3, 2), new RGB(10, 6, 4), new RGB(15, 9, 6), new RGB(20, 13, 8), new RGB(25, 16, 10), new RGB(30, 19, 12), new RGB(35, 22, 14), new RGB(40, 25, 16), new RGB(46, 28, 18), new RGB(51, 32, 20), new RGB(56, 35, 22), new RGB(61, 38, 24), new RGB(66, 41, 26), new RGB(71, 44, 28), new RGB(76, 47, 30), new RGB(81, 51, 32), new RGB(86, 54, 34), new RGB(91, 57, 36), new RGB(96, 60, 38), new RGB(101, 63, 40), new RGB(106, 66, 42), new RGB(111, 70, 44), new RGB(116, 73, 46), new RGB(121, 76, 48), new RGB(126, 79, 50), new RGB(132, 82, 52), new RGB(137, 85, 54), new RGB(142, 89, 56), new RGB(147, 92, 58), new RGB(152, 95, 60), new RGB(157, 98, 62), new RGB(162, 101, 64), new RGB(167, 104, 66), new RGB(172, 108, 68), new RGB(177, 111, 70), new RGB(182, 114, 72), new RGB(187, 117, 75), new RGB(192, 120, 77), new RGB(197, 123, 79), new RGB(202, 126, 81), new RGB(207, 130, 83), new RGB(212, 133, 85), new RGB(218, 136, 87), new RGB(223, 139, 89), new RGB(228, 142, 91), new RGB(233, 145, 93), new RGB(238, 149, 95), new RGB(243, 152, 97), new RGB(248, 155, 99), new RGB(253, 158, 101), new RGB(255, 161, 103), new RGB(255, 164, 105), new RGB(255, 168, 107), new RGB(255, 171, 109), new RGB(255, 174, 111), new RGB(255, 177, 113), new RGB(255, 180, 115), new RGB(255, 183, 117), new RGB(255, 187, 119), new RGB(255, 190, 121), new RGB(255, 193, 123), new RGB(255, 196, 125), new RGB(255, 199, 127) };
            }
        }
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] HSB
        {
            get
            {
                return new RGB[] { new RGB(255, 0, 0), new RGB(255, 24, 0), new RGB(255, 48, 0), new RGB(255, 72, 0), new RGB(255, 96, 0), new RGB(255, 120, 0), new RGB(255, 143, 0), new RGB(255, 167, 0), new RGB(255, 191, 0), new RGB(255, 215, 0), new RGB(255, 239, 0), new RGB(247, 255, 0), new RGB(223, 255, 0), new RGB(199, 255, 0), new RGB(175, 255, 0), new RGB(151, 255, 0), new RGB(128, 255, 0), new RGB(104, 255, 0), new RGB(80, 255, 0), new RGB(56, 255, 0), new RGB(32, 255, 0), new RGB(8, 255, 0), new RGB(0, 255, 16), new RGB(0, 255, 40), new RGB(0, 255, 64), new RGB(0, 255, 88), new RGB(0, 255, 112), new RGB(0, 255, 135), new RGB(0, 255, 159), new RGB(0, 255, 183), new RGB(0, 255, 207), new RGB(0, 255, 231), new RGB(0, 255, 255), new RGB(0, 231, 255), new RGB(0, 207, 255), new RGB(0, 183, 255), new RGB(0, 159, 255), new RGB(0, 135, 255), new RGB(0, 112, 255), new RGB(0, 88, 255), new RGB(0, 64, 255), new RGB(0, 40, 255), new RGB(0, 16, 255), new RGB(8, 0, 255), new RGB(32, 0, 255), new RGB(56, 0, 255), new RGB(80, 0, 255), new RGB(104, 0, 255), new RGB(128, 0, 255), new RGB(151, 0, 255), new RGB(175, 0, 255), new RGB(199, 0, 255), new RGB(223, 0, 255), new RGB(247, 0, 255), new RGB(255, 0, 239), new RGB(255, 0, 215), new RGB(255, 0, 191), new RGB(255, 0, 167), new RGB(255, 0, 143), new RGB(255, 0, 120), new RGB(255, 0, 96), new RGB(255, 0, 72), new RGB(255, 0, 48), new RGB(255, 0, 24) };
            }
        }
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] Jet
        {
            get
            {
                return new RGB[] { new RGB(0, 0, 143), new RGB(0, 0, 159), new RGB(0, 0, 175), new RGB(0, 0, 191), new RGB(0, 0, 207), new RGB(0, 0, 223), new RGB(0, 0, 239), new RGB(0, 0, 255), new RGB(0, 16, 255), new RGB(0, 32, 255), new RGB(0, 48, 255), new RGB(0, 64, 255), new RGB(0, 80, 255), new RGB(0, 96, 255), new RGB(0, 112, 255), new RGB(0, 128, 255), new RGB(0, 143, 255), new RGB(0, 159, 255), new RGB(0, 175, 255), new RGB(0, 191, 255), new RGB(0, 207, 255), new RGB(0, 223, 255), new RGB(0, 239, 255), new RGB(0, 255, 255), new RGB(16, 255, 239), new RGB(32, 255, 223), new RGB(48, 255, 207), new RGB(64, 255, 191), new RGB(80, 255, 175), new RGB(96, 255, 159), new RGB(112, 255, 143), new RGB(128, 255, 128), new RGB(143, 255, 112), new RGB(159, 255, 96), new RGB(175, 255, 80), new RGB(191, 255, 64), new RGB(207, 255, 48), new RGB(223, 255, 32), new RGB(239, 255, 16), new RGB(255, 255, 0), new RGB(255, 239, 0), new RGB(255, 223, 0), new RGB(255, 207, 0), new RGB(255, 191, 0), new RGB(255, 175, 0), new RGB(255, 159, 0), new RGB(255, 143, 0), new RGB(255, 128, 0), new RGB(255, 112, 0), new RGB(255, 96, 0), new RGB(255, 80, 0), new RGB(255, 64, 0), new RGB(255, 48, 0), new RGB(255, 32, 0), new RGB(255, 16, 0), new RGB(255, 0, 0), new RGB(239, 0, 0), new RGB(223, 0, 0), new RGB(207, 0, 0), new RGB(191, 0, 0), new RGB(175, 0, 0), new RGB(159, 0, 0), new RGB(143, 0, 0), new RGB(128, 0, 0) };
            }
        }
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] Pink
        {
            get
            {
                return new RGB[] { new RGB(30, 0, 0), new RGB(50, 26, 26), new RGB(64, 37, 37), new RGB(75, 45, 45), new RGB(85, 52, 52), new RGB(94, 59, 59), new RGB(102, 64, 64), new RGB(110, 69, 69), new RGB(117, 74, 74), new RGB(123, 79, 79), new RGB(130, 83, 83), new RGB(136, 87, 87), new RGB(141, 91, 91), new RGB(147, 95, 95), new RGB(152, 98, 98), new RGB(157, 102, 102), new RGB(162, 105, 105), new RGB(167, 108, 108), new RGB(172, 111, 111), new RGB(176, 114, 114), new RGB(181, 117, 117), new RGB(185, 120, 120), new RGB(189, 123, 123), new RGB(194, 126, 126), new RGB(195, 132, 129), new RGB(197, 138, 131), new RGB(199, 144, 134), new RGB(201, 149, 136), new RGB(202, 154, 139), new RGB(204, 159, 141), new RGB(206, 164, 144), new RGB(207, 169, 146), new RGB(209, 174, 148), new RGB(211, 178, 151), new RGB(212, 183, 153), new RGB(214, 187, 155), new RGB(216, 191, 157), new RGB(217, 195, 160), new RGB(219, 199, 162), new RGB(220, 203, 164), new RGB(222, 207, 166), new RGB(223, 211, 168), new RGB(225, 215, 170), new RGB(226, 218, 172), new RGB(228, 222, 174), new RGB(229, 225, 176), new RGB(231, 229, 178), new RGB(232, 232, 180), new RGB(234, 234, 185), new RGB(235, 235, 191), new RGB(237, 237, 196), new RGB(238, 238, 201), new RGB(240, 240, 206), new RGB(241, 241, 211), new RGB(243, 243, 216), new RGB(244, 244, 221), new RGB(245, 245, 225), new RGB(247, 247, 230), new RGB(248, 248, 234), new RGB(250, 250, 238), new RGB(251, 251, 243), new RGB(252, 252, 247), new RGB(254, 254, 251), new RGB(255, 255, 255) };
            }
        }
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] Autumn
        {
            get
            {
                return new RGB[] { new RGB(255, 0, 0), new RGB(255, 4, 0), new RGB(255, 8, 0), new RGB(255, 12, 0), new RGB(255, 16, 0), new RGB(255, 20, 0), new RGB(255, 24, 0), new RGB(255, 28, 0), new RGB(255, 32, 0), new RGB(255, 36, 0), new RGB(255, 40, 0), new RGB(255, 45, 0), new RGB(255, 49, 0), new RGB(255, 53, 0), new RGB(255, 57, 0), new RGB(255, 61, 0), new RGB(255, 65, 0), new RGB(255, 69, 0), new RGB(255, 73, 0), new RGB(255, 77, 0), new RGB(255, 81, 0), new RGB(255, 85, 0), new RGB(255, 89, 0), new RGB(255, 93, 0), new RGB(255, 97, 0), new RGB(255, 101, 0), new RGB(255, 105, 0), new RGB(255, 109, 0), new RGB(255, 113, 0), new RGB(255, 117, 0), new RGB(255, 121, 0), new RGB(255, 125, 0), new RGB(255, 130, 0), new RGB(255, 134, 0), new RGB(255, 138, 0), new RGB(255, 142, 0), new RGB(255, 146, 0), new RGB(255, 150, 0), new RGB(255, 154, 0), new RGB(255, 158, 0), new RGB(255, 162, 0), new RGB(255, 166, 0), new RGB(255, 170, 0), new RGB(255, 174, 0), new RGB(255, 178, 0), new RGB(255, 182, 0), new RGB(255, 186, 0), new RGB(255, 190, 0), new RGB(255, 194, 0), new RGB(255, 198, 0), new RGB(255, 202, 0), new RGB(255, 206, 0), new RGB(255, 210, 0), new RGB(255, 215, 0), new RGB(255, 219, 0), new RGB(255, 223, 0), new RGB(255, 227, 0), new RGB(255, 231, 0), new RGB(255, 235, 0), new RGB(255, 239, 0), new RGB(255, 243, 0), new RGB(255, 247, 0), new RGB(255, 251, 0), new RGB(255, 255, 0) };
            }
        }
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] Spring
        {
            get
            {
                return new RGB[] { new RGB(255, 0, 255), new RGB(255, 4, 251), new RGB(255, 8, 247), new RGB(255, 12, 243), new RGB(255, 16, 239), new RGB(255, 20, 235), new RGB(255, 24, 231), new RGB(255, 28, 227), new RGB(255, 32, 223), new RGB(255, 36, 219), new RGB(255, 40, 215), new RGB(255, 45, 210), new RGB(255, 49, 206), new RGB(255, 53, 202), new RGB(255, 57, 198), new RGB(255, 61, 194), new RGB(255, 65, 190), new RGB(255, 69, 186), new RGB(255, 73, 182), new RGB(255, 77, 178), new RGB(255, 81, 174), new RGB(255, 85, 170), new RGB(255, 89, 166), new RGB(255, 93, 162), new RGB(255, 97, 158), new RGB(255, 101, 154), new RGB(255, 105, 150), new RGB(255, 109, 146), new RGB(255, 113, 142), new RGB(255, 117, 138), new RGB(255, 121, 134), new RGB(255, 125, 130), new RGB(255, 130, 125), new RGB(255, 134, 121), new RGB(255, 138, 117), new RGB(255, 142, 113), new RGB(255, 146, 109), new RGB(255, 150, 105), new RGB(255, 154, 101), new RGB(255, 158, 97), new RGB(255, 162, 93), new RGB(255, 166, 89), new RGB(255, 170, 85), new RGB(255, 174, 81), new RGB(255, 178, 77), new RGB(255, 182, 73), new RGB(255, 186, 69), new RGB(255, 190, 65), new RGB(255, 194, 61), new RGB(255, 198, 57), new RGB(255, 202, 53), new RGB(255, 206, 49), new RGB(255, 210, 45), new RGB(255, 215, 40), new RGB(255, 219, 36), new RGB(255, 223, 32), new RGB(255, 227, 28), new RGB(255, 231, 24), new RGB(255, 235, 20), new RGB(255, 239, 16), new RGB(255, 243, 12), new RGB(255, 247, 8), new RGB(255, 251, 4), new RGB(255, 255, 0) };
            }
        }
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] Summer
        {
            get
            {
                return new RGB[] { new RGB(0, 128, 102), new RGB(4, 130, 102), new RGB(8, 132, 102), new RGB(12, 134, 102), new RGB(16, 136, 102), new RGB(20, 138, 102), new RGB(24, 140, 102), new RGB(28, 142, 102), new RGB(32, 144, 102), new RGB(36, 146, 102), new RGB(40, 148, 102), new RGB(45, 150, 102), new RGB(49, 152, 102), new RGB(53, 154, 102), new RGB(57, 156, 102), new RGB(61, 158, 102), new RGB(65, 160, 102), new RGB(69, 162, 102), new RGB(73, 164, 102), new RGB(77, 166, 102), new RGB(81, 168, 102), new RGB(85, 170, 102), new RGB(89, 172, 102), new RGB(93, 174, 102), new RGB(97, 176, 102), new RGB(101, 178, 102), new RGB(105, 180, 102), new RGB(109, 182, 102), new RGB(113, 184, 102), new RGB(117, 186, 102), new RGB(121, 188, 102), new RGB(125, 190, 102), new RGB(130, 192, 102), new RGB(134, 194, 102), new RGB(138, 196, 102), new RGB(142, 198, 102), new RGB(146, 200, 102), new RGB(150, 202, 102), new RGB(154, 204, 102), new RGB(158, 206, 102), new RGB(162, 208, 102), new RGB(166, 210, 102), new RGB(170, 212, 102), new RGB(174, 215, 102), new RGB(178, 217, 102), new RGB(182, 219, 102), new RGB(186, 221, 102), new RGB(190, 223, 102), new RGB(194, 225, 102), new RGB(198, 227, 102), new RGB(202, 229, 102), new RGB(206, 231, 102), new RGB(210, 233, 102), new RGB(215, 235, 102), new RGB(219, 237, 102), new RGB(223, 239, 102), new RGB(227, 241, 102), new RGB(231, 243, 102), new RGB(235, 245, 102), new RGB(239, 247, 102), new RGB(243, 249, 102), new RGB(247, 251, 102), new RGB(251, 253, 102), new RGB(255, 255, 102) };
            }
        }
        /// <summary>
        /// Returns the color scheme.
        /// </summary>
        public static RGB[] Winter
        {
            get
            {
                return new RGB[] { new RGB(0, 0, 255), new RGB(0, 4, 253), new RGB(0, 8, 251), new RGB(0, 12, 249), new RGB(0, 16, 247), new RGB(0, 20, 245), new RGB(0, 24, 243), new RGB(0, 28, 241), new RGB(0, 32, 239), new RGB(0, 36, 237), new RGB(0, 40, 235), new RGB(0, 45, 233), new RGB(0, 49, 231), new RGB(0, 53, 229), new RGB(0, 57, 227), new RGB(0, 61, 225), new RGB(0, 65, 223), new RGB(0, 69, 221), new RGB(0, 73, 219), new RGB(0, 77, 217), new RGB(0, 81, 215), new RGB(0, 85, 213), new RGB(0, 89, 210), new RGB(0, 93, 208), new RGB(0, 97, 206), new RGB(0, 101, 204), new RGB(0, 105, 202), new RGB(0, 109, 200), new RGB(0, 113, 198), new RGB(0, 117, 196), new RGB(0, 121, 194), new RGB(0, 125, 192), new RGB(0, 130, 190), new RGB(0, 134, 188), new RGB(0, 138, 186), new RGB(0, 142, 184), new RGB(0, 146, 182), new RGB(0, 150, 180), new RGB(0, 154, 178), new RGB(0, 158, 176), new RGB(0, 162, 174), new RGB(0, 166, 172), new RGB(0, 170, 170), new RGB(0, 174, 168), new RGB(0, 178, 166), new RGB(0, 182, 164), new RGB(0, 186, 162), new RGB(0, 190, 160), new RGB(0, 194, 158), new RGB(0, 198, 156), new RGB(0, 202, 154), new RGB(0, 206, 152), new RGB(0, 210, 150), new RGB(0, 215, 148), new RGB(0, 219, 146), new RGB(0, 223, 144), new RGB(0, 227, 142), new RGB(0, 231, 140), new RGB(0, 235, 138), new RGB(0, 239, 136), new RGB(0, 243, 134), new RGB(0, 247, 132), new RGB(0, 251, 130), new RGB(0, 255, 128) };
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model RYB.
    /// </summary>
    [Serializable]
    public struct RYB : IColorSpace, ICloneable
    {
        #region Private data
        private byte r;
        private byte y;
        private byte b;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure RYB.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="yellow">Yellow [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        public RYB(int red, int yellow, int blue)
        {
            this.r = (byte)((red > 255) ? 255 : ((red < 0) ? 0 : red));
            this.y = (byte)((yellow > 255) ? 255 : ((yellow < 0) ? 0 : yellow));
            this.b = (byte)((blue > 255) ? 255 : ((blue < 0) ? 0 : blue));
        }
        /// <summary>
        /// Creates an instance of the structure RYB.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="yellow">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        public RYB(double red, double yellow, double blue)
        {
            this.r = (byte)((red > 255) ? 255 : ((red < 0) ? 0 : red));
            this.y = (byte)((yellow > 255) ? 255 : ((yellow < 0) ? 0 : yellow));
            this.b = (byte)((blue > 255) ? 255 : ((blue < 0) ? 0 : blue));
        }
        /// <summary>
        /// Defines a component of the color model [0, 255].
        /// </summary>
        public byte Red
        {
            get
            {
                return r;
            }
            set
            {
                r = (byte)((value > 255) ? 255 : ((value < 0) ? 0 : value));
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 255].
        /// </summary>
        public byte Yellow
        {
            get
            {
                return y;
            }
            set
            {
                y = (byte)((value > 255) ? 255 : ((value < 0) ? 0 : value));
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 255].
        /// </summary>
        public byte Blue
        {
            get
            {
                return b;
            }
            set
            {
                b = (byte)((value > 255) ? 255 : ((value < 0) ? 0 : value));
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">RYB structure</param>
        /// <param name="item2">RYB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(RYB item1, RYB item2)
        {
            return (
                item1.Red == item2.Red
                && item1.Yellow == item2.Yellow
                && item1.Blue == item2.Blue
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">RYB structure</param>
        /// <param name="item2">RYB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(RYB item1, RYB item2)
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

            return (this == (RYB)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Red.GetHashCode() ^ Yellow.GetHashCode() ^ Blue.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Red.ToString() + "\n" + Yellow.ToString() + "\n" + Blue.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new RYB(this.Red, this.Yellow, this.Blue);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public RYB Clone()
        {
            return new RYB(this.Red, this.Yellow, this.Blue);
        }
        #endregion

        #region RYB convert
        /// <summary>
        /// Converts a color model RGB in model RYB.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>RYB structure</returns>
        public static RYB FromRGB(int red, int green, int blue)
        {
            // Arah J. Leonard
            // https://github.com/bahamas10/node-rgb2ryb/blob/master/rgb2ryb.js
            // http://www.deathbysoftware.com/colors/index.html
            //

            double r = red, g = green, b = blue;
            double w = Maths.Min(r, g, b);

            r -= w; g -= w; b -= w;

            double mg = Maths.Max(r, g, b);
            double y = Maths.Min(r, g);

            r -= y;
            g -= y;

            //if (b != g)
            {
                b /= 2.0;
                g /= 2.0;
            }

            y += g;
            b += g;

            double my = Maths.Max(r, y, b);

            //if (my > 0)
            {
                double n = mg / my;
                r *= n;
                y *= n;
                b *= n;
            }

            // Add the white back in.
            r += w;
            y += w;
            b += w;

            return new RYB(r, y, b);
        }
        /// <summary>
        /// Converts a color model RGB in model RYB.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>RYB structure</returns>
        public static RYB FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model RYB in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                // Arah J. Leonard
                // https://github.com/bahamas10/node-rgb2ryb/blob/master/rgb2ryb.js
                // http://www.deathbysoftware.com/colors/index.html
                //

                double rr = r;
                double yy = y;
                double bb = b;

                double w = Maths.Min(rr, yy, bb);
                rr -= w;
                yy -= w;
                bb -= w;

                double my = Maths.Max(rr, yy, bb);

                // Get the green out of the yellow and blue
                double gg = Maths.Min(yy, bb);
                yy -= gg;
                bb -= gg;

                bb *= 2.0;
                gg *= 2.0;

                // Redistribute the remaining yellow.
                rr += yy;
                gg += yy;

                // Normalize to values.
                double mg = Maths.Max(rr, gg, bb);
                double n = my / mg;
                rr *= n;
                gg *= n;
                bb *= n;

                // Add the white back in.
                rr += w;
                gg += w;
                bb += w;

                return new RGB(rr, gg, bb);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model sRGB.
    /// </summary>
    [Serializable]
    public struct sRGB : IColorSpace, ICloneable
    {
        #region private data
        private double r;
        private double g;
        private double b;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure sRGB.
        /// </summary>
        /// <param name="red">Red [0, 1]</param>
        /// <param name="green">Green [0, 1]</param>
        /// <param name="blue">Blue [0, 1]</param>
        public sRGB(double red, double green, double blue)
        {
            this.r = (red > 1) ? 1 : ((red < 0) ? 0 : red);
            this.g = (green > 1) ? 1 : ((green < 0) ? 0 : green);
            this.b = (blue > 1) ? 1 : ((blue < 0) ? 0 : blue);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Red
        {
            get
            {
                return r;
            }
            set
            {
                r = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Green
        {
            get
            {
                return g;
            }
            set
            {
                g = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Blue
        {
            get
            {
                return b;
            }
            set
            {
                b = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">sRGB structure</param>
        /// <param name="item2">sRGB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(sRGB item1, sRGB item2)
        {
            return (
                item1.Red == item2.Red
                && item1.Green == item2.Green
                && item1.Blue == item2.Blue
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">sRGB structure</param>
        /// <param name="item2">sRGB structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(sRGB item1, sRGB item2)
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

            return (this == (sRGB)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Red.GetHashCode() ^ Green.GetHashCode() ^ Blue.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Red.ToString() + "\n" + Green.ToString() + "\n" + Blue.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new sRGB(this.Red, this.Green, this.Blue);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public sRGB Clone()
        {
            return new sRGB(this.Red, this.Green, this.Blue);
        }
        #endregion

        #region sRGB convert
        /// <summary>
        /// Converts a color model RGB in model sRGB.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>sRGB structure</returns>
        public static sRGB FromRGB(int red, int green, int blue)
        {
            double r = red / 255.0f;
            double g = green / 255.0f;
            double b = blue / 255.0f;

            return new sRGB(r, g, b);
        }
        /// <summary>
        /// Converts a color model RGB in model sRGB.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>sRGB structure</returns>
        public static sRGB FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model sRGB in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                return new RGB(r * 255.0, g * 255.0, b * 255.0);

            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model YUV.
    /// </summary>
    [Serializable]
    public struct YUV : IColorSpace, ICloneable
    {
        #region Private data
        private double y;
        private double u;
        private double v;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YUV.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="u">U [-0.436, 0.436]</param>
        /// <param name="v">V [-0.614, 0.614]</param>
        public YUV(double y, double u, double v)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.u = (u > 0.436f) ? 0.436f : ((u < -0.436f) ? -0.436f : u);
            this.v = (v > 0.614f) ? 0.614f : ((v < -0.614f) ? -0.614f : v);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Y
        {
            get
            {
                return y;
            }
            set
            {
                y = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.436, 0.436].
        /// </summary>
        public double U
        {
            get
            {
                return u;
            }
            set
            {
                u = (value > 0.436f) ? 0.436f : ((value < -0.436f) ? -0.436f : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.614, 0.614].
        /// </summary>
        public double V
        {
            get
            {
                return v;
            }
            set
            {
                v = (value > 0.614f) ? 0.614f : ((value < -0.614f) ? -0.614f : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YUV structure</param>
        /// <param name="item2">YUV structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YUV item1, YUV item2)
        {
            return (
                item1.Y == item2.Y
                && item1.U == item2.U
                && item1.V == item2.V
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YUV structure</param>
        /// <param name="item2">YUV structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YUV item1, YUV item2)
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

            return (this == (YUV)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ U.GetHashCode() ^ V.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + U.ToString() + "\n" + V.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YUV(this.Y, this.U, this.V);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YUV Clone()
        {
            return new YUV(this.Y, this.U, this.V);
        }
        #endregion

        #region YUV convert
        /// <summary>
        /// Converts a color model RGB in model YUV.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YUV structure</returns>
        public static YUV FromRGB(int red, int green, int blue)
        {
            double r = red / 255.0f;
            double g = green / 255.0f;
            double b = blue / 255.0f;

            double Y = 0.299f * r + 0.587f * g + 0.114f * b;
            double U = -0.147f * r - 0.288f * g + 0.436f * b;
            double V = 0.615f * r - 0.514f * g - 0.100f * b;

            return new YUV(Y, U, V);
        }
        /// <summary>
        /// Converts a color model RGB in model YUV.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YUV structure</returns>
        public static YUV FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YUV in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 1.139f * v) * 255);
                int g = (int)((y - 0.394f * u - 0.580f * v) * 255);
                int b = (int)((y + 2.032f * u) * 255);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model YIQ.
    /// </summary>
    [Serializable]
    public struct YIQ : IColorSpace, ICloneable
    {
        #region Private data
        private double y;
        private double i;
        private double q;
        #endregion

        #region Stucture components
        /// <summary>
        /// Creates an instance of the structure YIQ.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="i">I [-0.5957, 0.5957]</param>
        /// <param name="q">Q [-0.5226, 0.5226]</param>
        public YIQ(double y, double i, double q)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.i = (i > 0.5957f) ? 0.5957f : ((i < -0.5957f) ? -0.5957f : i);
            this.q = (q > 0.5226f) ? 0.5226f : ((q < -0.5226f) ? -0.5226f : q);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Y
        {
            get
            {
                return y;
            }
            set
            {
                y = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.5957, 0.5957].
        /// </summary>
        public double I
        {
            get
            {
                return i;
            }
            set
            {
                i = (value > 0.5957f) ? 0.5957f : ((value < -0.5957f) ? -0.5957f : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.5226, 0.5226].
        /// </summary>
        public double Q
        {
            get
            {
                return q;
            }
            set
            {
                q = (value > 0.5226f) ? 0.5226f : ((value < -0.5226f) ? -0.5226f : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YIQ structure</param>
        /// <param name="item2">YIQ structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YIQ item1, YIQ item2)
        {
            return (
                item1.Y == item2.Y
                && item1.I == item2.I
                && item1.Q == item2.Q
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YIQ structure</param>
        /// <param name="item2">YIQ structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YIQ item1, YIQ item2)
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

            return (this == (YIQ)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ I.GetHashCode() ^ Q.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + I.ToString() + "\n" + Q.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YIQ(this.Y, this.I, this.Q);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YIQ Clone()
        {
            return new YIQ(this.Y, this.I, this.Q);
        }
        #endregion

        #region YIQ convert
        /// <summary>
        /// Converts a color model RGB in model YIQ.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YIQ structure</returns>
        public static YIQ FromRGB(int red, int green, int blue)
        {
            double r = red / 255.0f;
            double g = green / 255.0f;
            double b = blue / 255.0f;

            double Y = 0.299f * r + 0.587f * g + 0.114f * b;
            double I = 0.596f * r - 0.274f * g - 0.322f * b;
            double Q = 0.211f * r - 0.522f * g + 0.311f * b;

            return new YIQ(Y, I, Q);
        }
        /// <summary>
        /// Converts a color model RGB in model YIQ.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YIQ structure</returns>
        public static YIQ FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YIQ in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 0.956f * i + 0.623f * q) * 255);
                int g = (int)((y - 0.272f * i - 0.648f * q) * 255);
                int b = (int)((y - 1.105f * i + 1.705f * q) * 255);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model YCbCr.
    /// </summary>
    [Serializable]
    public struct YCbCr : IColorSpace, ICloneable
    {
        #region Private data
        private double y;
        private double cb;
        private double cr;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YCbCr.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="cb">Cb [-1, 1]</param>
        /// <param name="cr">Cr [-1, 1]</param>
        public YCbCr(double y, double cb, double cr)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.cb = (cb > 1) ? 1 : ((cb < -1) ? -1 : cb);
            this.cr = (cr > 1) ? 1 : ((cr < -1) ? -1 : cr);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Y
        {
            get
            {
                return y;
            }
            set
            {
                y = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-1, 1].
        /// </summary>
        public double Cb
        {
            get
            {
                return cb;
            }
            set
            {
                cb = (value > 1) ? 1 : ((value < -1) ? -1 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-1, 1].
        /// </summary>
        public double Cr
        {
            get
            {
                return cr;
            }
            set
            {
                cr = (value > 1) ? 1 : ((value < -1) ? -1 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YCbCr structure</param>
        /// <param name="item2">YCbCr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YCbCr item1, YCbCr item2)
        {
            return (
                item1.Y == item2.Y
                && item1.Cb == item2.Cb
                && item1.Cr == item2.Cr
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YCbCr structure</param>
        /// <param name="item2">YCbCr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YCbCr item1, YCbCr item2)
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

            return (this == (YCbCr)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ Cb.GetHashCode() ^ Cr.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + Cb.ToString() + "\n" + Cr.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YCbCr(this.Y, this.Cb, this.Cr);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YCbCr Clone()
        {
            return new YCbCr(this.Y, this.Cb, this.Cr);
        }
        #endregion

        #region YCbCr convert
        /// <summary>
        /// Converts a color model RGB in model YCbCr.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YCbCr structure</returns>
        public static YCbCr FromRGB(int red, int green, int blue)
        {
            double r = red / 255.0;
            double g = green / 255.0;
            double b = blue / 255.0;

            double Y = 0.299 * r + 0.587 * g + 0.114 * b;
            double Cb = -0.172 * r - 0.339 * g + 0.511 * b + 0.5;
            double Cr = 0.511 * r - 0.428 * g - 0.083 * b + 0.5;

            return new YCbCr(Y, Cb, Cr);
        }
        /// <summary>
        /// Converts a color model RGB in model YCbCr.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YCbCr structure</returns>
        public static YCbCr FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YCbCr in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 1.371f * (cr - 0.5f)) * 255);
                int g = (int)((y - 0.698f * (cr - 0.5f) - 0.336f * (cb - 0.5f)) * 255);
                int b = (int)((y + 1.732f * (cb - 0.5f)) * 255);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model YDbDr.
    /// </summary>
    [Serializable]
    public struct YDbDr : IColorSpace, ICloneable
    {
        #region Private data
        private double y;
        private double db;
        private double dr;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YDbDr.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="db">Db [-1.333, 1.333]</param>
        /// <param name="dr">Dr [-1.333, 1.333]</param>
        public YDbDr(double y, double db, double dr)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.db = (db > 1.333) ? 1.333 : ((db < -1.333) ? -1.333 : db);
            this.dr = (dr > 1.333) ? 1.333 : ((dr < -1.333) ? -1.333 : dr);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Y
        {
            get
            {
                return y;
            }
            set
            {
                y = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-1.333, 1.333].
        /// </summary>
        public double Db
        {
            get
            {
                return db;
            }
            set
            {
                db = (value > 1.333) ? 1 : ((value < -1.333) ? -1.333 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-1.333, 1.333].
        /// </summary>
        public double Dr
        {
            get
            {
                return dr;
            }
            set
            {
                dr = (value > 1.333) ? 1.333 : ((value < -1.333) ? -1.333 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YDbDr structure</param>
        /// <param name="item2">YDbDr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YDbDr item1, YDbDr item2)
        {
            return (
                item1.Y == item2.Y
                && item1.Db == item2.Db
                && item1.Dr == item2.Dr
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YDbDr structure</param>
        /// <param name="item2">YDbDr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YDbDr item1, YDbDr item2)
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

            return (this == (YDbDr)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ Db.GetHashCode() ^ Dr.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + Db.ToString() + "\n" + Dr.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YDbDr(this.Y, this.Db, this.Dr);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YDbDr Clone()
        {
            return new YDbDr(this.Y, this.Db, this.Dr);
        }
        #endregion

        #region YDbDr convert
        /// <summary>
        /// Converts a color model RGB in model YDbDr.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YDbDr structure</returns>
        public static YDbDr FromRGB(int red, int green, int blue)
        {
            double r = red / 255.0f;
            double g = green / 255.0f;
            double b = blue / 255.0f;

            double Y = 0.299f * r + 0.587f * g + 0.114f * b;
            double Cb = -0.450f * r - 0.883f * g + 1.333f * b;
            double Cr = -1.333f * r - 1.116f * g - 0.217f * b;

            return new YDbDr(Y, Cb, Cr);
        }
        /// <summary>
        /// Converts a color model RGB in model YDbDr.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YDbDr structure</returns>
        public static YDbDr FromRGB(RGB rgb)
        {
            return YDbDr.FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YDbDr in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 0.00009 * db - 0.52591 * dr) * 255);
                int g = (int)((y - 0.12913 * db + 0.26789 * dr) * 255);
                int b = (int)((y + 0.66467 * db - 0.00007 * dr) * 255);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model YCgCo.
    /// </summary>
    [Serializable]
    public struct YCgCo : IColorSpace, ICloneable
    {
        #region Private data
        private double y;
        private double cg;
        private double co;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YDbDr.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="cg">Cg [-0.5, 0.5]</param>
        /// <param name="co">Co [-0.5, 0.5]</param>
        public YCgCo(double y, double cg, double co)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.cg = (cg > 0.5) ? 0.5 : ((cg < -0.5) ? -0.5 : cg);
            this.co = (co > 0.5) ? 0.5 : ((co < -0.5) ? -0.5 : co);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Y
        {
            get
            {
                return y;
            }
            set
            {
                y = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.5, 0.5].
        /// </summary>
        public double Cg
        {
            get
            {
                return cg;
            }
            set
            {
                cg = (value > 0.5) ? 0.5: ((value < -0.5) ? -0.5 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.5, 0.5].
        /// </summary>
        public double Co
        {
            get
            {
                return co;
            }
            set
            {
                co = (value > 0.5) ? 0.5 : ((value < -0.5) ? -0.5 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YCgCo structure</param>
        /// <param name="item2">YCgCo structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YCgCo item1, YCgCo item2)
        {
            return (
                item1.Y == item2.Y
                && item1.Cg == item2.Cg
                && item1.Co == item2.Co
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YCgCo structure</param>
        /// <param name="item2">YCgCo structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YCgCo item1, YCgCo item2)
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

            return (this == (YCgCo)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ Cg.GetHashCode() ^ Co.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + Cg.ToString() + "\n" + Co.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YCgCo(this.Y, this.Cg, this.Co);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YCgCo Clone()
        {
            return new YCgCo(this.Y, this.Cg, this.Co);
        }
        #endregion

        #region YCgCo convert
        /// <summary>
        /// Converts a color model RGB in model YCgCo.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YCgCo structure</returns>
        public static YCgCo FromRGB(int red, int green, int blue)
        {
            double r = red / 255.0;
            double g = green / 255.0;
            double b = blue / 255.0;

            double Y = 0.25 * r + 0.5 * g + 0.25 * b;
            double Cg = -0.25 * r + 0.5 * g - 0.25 * b;
            double Co = 0.5 * r - 0.0 * g - 0.5 * b;

            return new YCgCo(Y, Cg, Co);
        }
        /// <summary>
        /// Converts a color model RGB in model YCgCo.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YCgCo structure</returns>
        public static YCgCo FromRGB(RGB rgb)
        {
            return FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YCgCo in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y - cg + co) * 255.0);
                int g = (int)((y + cg) * 255.0);
                int b = (int)((y - cg - co) * 255.0);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines a color model YPbPr.
    /// </summary>
    [Serializable]
    public struct YPbPr : IColorSpace, ICloneable
    {
        #region Private data
        private double y;
        private double pb;
        private double pr;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure YPbPr.
        /// </summary>
        /// <param name="y">Y [0, 1]</param>
        /// <param name="pb">Pb [-0.5, 0.5]</param>
        /// <param name="pr">Pr [-0.5, 0.5]</param>
        public YPbPr(double y, double pb, double pr)
        {
            this.y = (y > 1) ? 1 : ((y < 0) ? 0 : y);
            this.pb = (pb > 0.5) ? 0.5 : ((pb < -0.5) ? -0.5 : pb);
            this.pr = (pr > 0.5) ? 0.5 : ((pr < -0.5) ? -0.5 : pr);
        }
        /// <summary>
        /// Defines a component of the color model [0, 1].
        /// </summary>
        public double Y
        {
            get
            {
                return y;
            }
            set
            {
                y = (value > 1) ? 1 : ((value < 0) ? 0 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.5, 0.5].
        /// </summary>
        public double Pb
        {
            get
            {
                return pb;
            }
            set
            {
                pb = (value > 0.5) ? 0.5 : ((value < -0.5) ? -0.5 : value);
            }
        }
        /// <summary>
        /// Defines a component of the color model [-0.5, 0.5].
        /// </summary>
        public double Pr
        {
            get
            {
                return pr;
            }
            set
            {
                pr = (value > 0.5) ? 0.5 : ((value < -0.5) ? -0.5 : value);
            }
        }
        #endregion

        #region Boolean
        /// <summary>
        /// Checks the equality of two class objects.
        /// </summary>
        /// <param name="item1">YPbPr structure</param>
        /// <param name="item2">YPbPr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator ==(YPbPr item1, YPbPr item2)
        {
            return (
                item1.Y == item2.Y
                && item1.Pb == item2.Pb
                && item1.Pr == item2.Pr
                );
        }
        /// <summary>
        /// Checks the inequality of two class objects.
        /// </summary>
        /// <param name="item1">YPbPr structure</param>
        /// <param name="item2">YPbPr structure</param>
        /// <returns>Boolean</returns>
        public static bool operator !=(YPbPr item1, YPbPr item2)
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

            return (this == (YPbPr)obj);
        }
        /// <summary>
        /// Plays the role of a hash function of a certain type.
        /// </summary>
        /// <returns>Integer number</returns>
        public override int GetHashCode()
        {
            return Y.GetHashCode() ^ Pb.GetHashCode() ^ Pr.GetHashCode();
        }
        /// <summary>
        /// Returns a System.String object that represents the current object.
        /// </summary>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public override string ToString()
        {
            return Y.ToString() + "\n" + Pb.ToString() + "\n" + Pr.ToString();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        object ICloneable.Clone()
        {
            return new YPbPr(this.Y, this.Pb, this.Pr);
        }
        /// <summary>
        /// Creates a copy of the color model.
        /// </summary>
        /// <returns>Structure</returns>
        public YPbPr Clone()
        {
            return new YPbPr(this.Y, this.Pb, this.Pr);
        }
        #endregion

        #region YPbPr convert
        /// <summary>
        /// Converts a color model RGB in model YPbPr.
        /// </summary>
        /// <param name="red">Red [0, 255]</param>
        /// <param name="green">Green [0, 255]</param>
        /// <param name="blue">Blue [0, 255]</param>
        /// <returns>YPbPr structure</returns>
        public static YPbPr FromRGB(int red, int green, int blue)
        {
            double r = red / 255.0;
            double g = green / 255.0;
            double b = blue / 255.0;

            double Y = 0.299 * r + 0.587 * g + 0.114 * b;
            double Cb = -0.169 * r - 0.331 * g + 0.500 * b;
            double Cr = 0.500 * r - 0.419 * g - 0.081 * b;

            return new YPbPr(Y, Cb, Cr);
        }
        /// <summary>
        /// Converts a color model RGB in model YPbPr.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>YPbPr structure</returns>
        public static YPbPr FromRGB(RGB rgb)
        {
            return YPbPr.FromRGB(rgb.Red, rgb.Green, rgb.Blue);
        }
        #endregion

        #region RGB convert
        /// <summary>
        /// Converts a color model YPbPr in model RGB.
        /// </summary>
        /// <returns>RGB structure</returns>
        public RGB ToRGB
        {
            get
            {
                int r = (int)((y + 1.402 * pr) * 255.0);
                int g = (int)((y - 0.344 * pb - 0.714 * pr) * 255.0);
                int b = (int)((y + 1.772 * pb) * 255.0);

                return new RGB(r, g, b);
            }
        }
        #endregion
    }
    #endregion

    #region Unknown colorspace
    /// <summary>
    /// Defines an unknown color model.
    /// This color model can play the role of any color space.
    /// </summary>
	[Serializable]
    public struct Unknown : IColorSpace, ICloneable
    {
        #region Private data
        private double x;
        private double y;
        private double z;
        #endregion

        #region Structure components
        /// <summary>
        /// Creates an instance of the structure.
        /// </summary>
        /// <param name="x">Component X</param>
        /// <param name="y">Component Y</param>
        /// <param name="z">Component Z</param>
        public Unknown(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        /// <summary>
        /// Defines the component of the color model.
        /// </summary>
        public double X
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
        public double Y
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
        public double Z
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
            return X.ToString() + "\n" + Y.ToString() + "\n" + Z.ToString();
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
    #endregion

    #region Interfaces
    /// <summary>
    /// Defines the color space interface.
    /// </summary>
    public interface IColorSpace
    {
        #region Interface
        /// <summary>
        /// Returns the color model RGB.
        /// </summary>
        RGB ToRGB { get; }
        #endregion
    }
    #endregion
}
