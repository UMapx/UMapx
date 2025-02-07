using System;
using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Colorspace
{
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
        public RGB(float red, float green, float blue)
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
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public static float Average(float red, float green, float blue)
        {
            return (red + green + blue) / 3.0f;
        }
        /// <summary>
        /// Calculates the average brightness value.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>float precision floating point number</returns>
        public static int Average(RGB rgb)
        {
            return Average(rgb.Red, rgb.Green, rgb.Blue);
        }
        /// <summary>
        /// Calculates the average brightness value.
        /// </summary>
        /// <param name="rgb">sRGB structure</param>
        /// <returns>float precision floating point number</returns>
        public static float Average(sRGB rgb)
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
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public static float PAL(float red, float green, float blue)
        {
            return 0.299f * red + 0.587f * green + 0.114f * blue;
        }
        /// <summary>
        /// Calculates the brightness value in the standard (PAL/NTC).
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>float precision floating point number</returns>
        public static float PAL(RGB rgb)
        {
            return PAL(rgb.Red, rgb.Green, rgb.Blue);
        }
        /// <summary>
        /// Calculates the brightness value in the standard (PAL/NTC).
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>float precision floating point number</returns>
        public static float PAL(sRGB rgb)
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
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public static float HDTV(float red, float green, float blue)
        {
            return 0.2126f * red + 0.7152f * green + 0.0722f * blue;
        }
        /// <summary>
        /// Calculates the brightness value in the standard HDTV.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>float precision floating point number</returns>
        public static int HDTV(RGB rgb)
        {
            return HDTV(rgb.Red, rgb.Green, rgb.Blue);
        }
        /// <summary>
        /// Calculates the brightness value in the standard HDTV.
        /// </summary>
        /// <param name="rgb">sRGB structure</param>
        /// <returns>float precision floating point number</returns>
        public static float HDTV(sRGB rgb)
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
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public static float RYY(float red, float green, float blue)
        {
            return 0.5f * red + 0.419f * green + 0.081f * blue;
        }
        /// <summary>
        /// Calculates the brightness value in the standard RYY.
        /// </summary>
        /// <param name="rgb">RGB structure</param>
        /// <returns>float precision floating point number</returns>
        public static int RYY(RGB rgb)
        {
            return RYY(rgb.Red, rgb.Green, rgb.Blue);
        }
        /// <summary>
        /// Calculates the brightness value in the standard RYY.
        /// </summary>
        /// <param name="rgb">sRGB structure</param>
        /// <returns>float precision floating point number</returns>
        public static float RYY(sRGB rgb)
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
        public static RGB Temp2RGB(float temperature)
        {
            // Approximation of Planckian locus in RGB model
            // Designed by Valery Asiryan, 2018

            float r = 0, g, b = 0;           // color channels,
            float x = temperature / 1000.0f; // normalize temerature,
            float y = x * x, z = y * x;      // variables.

            // Approximate red channel:
            if (x < 5)
            {
                r = 255;
            }
            else if (x <= 7)
            {
                r = -1575f + 718.5f * x - 70.5f * y;
            }
            // ************************************************************************

            // Approximate green channel:
            g = -50.1667f + 99.9846f * x - 7.7844f * y - 0.1175f * z;
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
        public static RGB Saturation(int red, int green, int blue, float s)
        {
            float max = Maths.Max(red, green, blue);

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
        public static sRGB Saturation(float red, float green, float blue, float s)
        {
            float max = Maths.Max(red, green, blue);

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
        public static RGB Saturation(RGB rgb, float s)
        {
            return Saturation(rgb.Red, rgb.Green, rgb.Blue, s);
        }
        /// <summary>
        /// Corrects color saturation.
        /// </summary>
        /// <param name="rgb">sRGB structure</param>
        /// <param name="s">Saturation</param>
        /// <returns>RGB structure</returns>
        public static sRGB Saturation(sRGB rgb, float s)
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
        public static RGB[] SchemeFromHue(float hue, uint length)
        {
            RGB[] scheme = new RGB[length];
            hue %= 360.0f;

            for (int i = 0; i < scheme.Length; i++)
            {
                scheme[i] = (new HSL(hue, 1.0f, i / (scheme.Length - 1.0f))).ToRGB;
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
}
