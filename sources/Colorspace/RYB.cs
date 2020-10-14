using System;
using UMapx.Core;

namespace UMapx.Colorspace
{
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
}
