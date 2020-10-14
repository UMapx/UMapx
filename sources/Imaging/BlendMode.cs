using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Used to blending layers.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.pegtop.net/delphi/articles/blendmodes/index.htm
    /// </remarks>
    /// </summary>
    public static class BlendMode
    {
        #region Blend mode components
        /// <summary>
        /// Implements the averaging function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Average(double a, double b)
        {
            return (a + b) / 2.0;
        }
        /// <summary>
        /// Implements the screening function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Screen(double a, double b)
        {
            return 1 - (1 - a) * (1 - b);
        }
        /// <summary>
        /// Implements the difference function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Difference(double a, double b)
        {
            return Math.Abs(a - b);
        }
        /// <summary>
        /// Implements the negation function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Negation(double a, double b)
        {
            return 1 - Math.Abs(1 - a - b);
        }
        /// <summary>
        /// Implements the exclusion function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Exclusion(double a, double b)
        {
            return a + b - 2 * a * b;
        }
        /// <summary>
        /// Implements the overlaying function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Overlay(double a, double b)
        {
            if (a < 0.5)
            {
                return 2 * a * b;
            }
            return 1 - 2 * (1 - a) * (1 - b);
        }
        /// <summary>
        /// Implements the "hard light" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double HardLight(double a, double b)
        {
            if (b < 0.5)
            {
                return 2 * a * b;
            }
            return 1 - 2 * (1 - a) * (1 - b);
        }
        /// <summary>
        /// Implements the "dodge" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Dodge(double a, double b)
        {
            return a / (1 - b);
        }
        /// <summary>
        /// Implements the "soft dodge" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double SoftDodge(double a, double b)
        {
            if (a + b < 1)
            {
                return 0.5 * a / (1 - b);
            }
            return 1 - 0.5 * (1 - b) / a;
        }
        /// <summary>
        /// Implements the "burn" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Burn(double a, double b)
        {
            return 1 - (1 - a) / b;
        }
        /// <summary>
        /// Implements the "soft burn" function".
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double SoftBurn(double a, double b)
        {
            if (a + b < 1)
            {
                return 0.5 * b / (1 - a);
            }
            return 1 - 0.5 * (1 - a) / b;
        }
        /// <summary>
        /// Implements the reflection function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Reflect(double a, double b)
        {
            return a * a / (1 - b);
        }
        /// <summary>
        /// Implements the glow function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Glow(double a, double b)
        {
            return b * b / (1 - a);
        }
        /// <summary>
        /// Implements the stamp function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Stamp(double a, double b)
        {
            return a + 2 * b - 1;
        }
        /// <summary>
        /// Implements the "freeze" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Freeze(double a, double b)
        {
            double x = 1 - a;
            return 1 - x * x / b;
        }
        /// <summary>
        /// Implements the "heat" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Heat(double a, double b)
        {
            double x = 1 - b;
            return 1 - x * x / a;
        }
        /// <summary>
        /// Implements the interpolation function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Interpolation(double a, double b)
        {
            return 0.5 - 0.25 * Math.Cos(Math.PI * a) - 0.25 * Math.Cos(Math.PI * b);
        }
        /// <summary>
        /// Implements the function of "soft light" (Adobe Photoshop).
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Photoshop(double a, double b)
        {
            if (b < 0.5)
            {
                return 2 * a * b + a * a * (1 - 2 * a);
            }
            return 2 * a * (1 - b) + Math.Sqrt(a) * (2 * b - 1);
        }
        /// <summary>
        /// Implements the function of "soft light" (Illusions.hu).
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Illusions(double a, double b)
        {
            double x = 2 * (0.5 - b);
            double y = Math.Pow(2, x);
            return Math.Pow(a, y);
        }
        /// <summary>
        /// Implements the function of "soft light" (Pegtop).
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Pegtop(double a, double b)
        {
            return (1 - 2 * b) * a * a + 2 * b * a;
        }
        /// <summary>
        /// Implements the "Cairo" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Fw3c(double a, double b)
        {
            if (b <= 0.5)
            {
                return a - (1 - 2 * b) * a * (1 - a);
            }
            return a + (2 * b - 1) * (BlendMode.Gw3c(a) - a);
        }
        /// <summary>
        /// Implements the "Cairo" function.
        /// </summary>
        /// <param name="a">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gw3c(double a)
        {
            if (a <= 0.25)
            {
                return ((16 * a - 12) * a + 4) * a;
            }
            return Math.Sqrt(a);
        }
        #endregion
    }
}
