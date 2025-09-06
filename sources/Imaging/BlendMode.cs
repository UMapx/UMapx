using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Used to blending layers.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.pegtop.net/delphi/articles/blendmodes/index.htm
    /// </remarks>
    public static class BlendMode
    {
        #region Blend mode components
        /// <summary>
        /// Implements the averaging function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Average(float a, float b)
        {
            return (a + b) / 2.0f;
        }
        /// <summary>
        /// Implements the screening function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Screen(float a, float b)
        {
            return 1 - (1 - a) * (1 - b);
        }
        /// <summary>
        /// Implements the difference function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Difference(float a, float b)
        {
            return Math.Abs(a - b);
        }
        /// <summary>
        /// Implements the negation function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Negation(float a, float b)
        {
            return 1 - Math.Abs(1 - a - b);
        }
        /// <summary>
        /// Implements the exclusion function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Exclusion(float a, float b)
        {
            return a + b - 2 * a * b;
        }
        /// <summary>
        /// Implements the overlaying function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Overlay(float a, float b)
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
        /// <returns>Value</returns>
        public static float HardLight(float a, float b)
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
        /// <returns>Value</returns>
        public static float Dodge(float a, float b)
        {
            return a / (1 - b);
        }
        /// <summary>
        /// Implements the "soft dodge" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float SoftDodge(float a, float b)
        {
            if (a + b < 1)
            {
                return 0.5f * a / (1 - b);
            }
            return 1 - 0.5f * (1 - b) / a;
        }
        /// <summary>
        /// Implements the "burn" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Burn(float a, float b)
        {
            return 1 - (1 - a) / b;
        }
        /// <summary>
        /// Implements the "soft burn" function".
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float SoftBurn(float a, float b)
        {
            if (a + b < 1)
            {
                return 0.5f * b / (1 - a);
            }
            return 1 - 0.5f * (1 - a) / b;
        }
        /// <summary>
        /// Implements the reflection function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Reflect(float a, float b)
        {
            return a * a / (1 - b);
        }
        /// <summary>
        /// Implements the glow function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Glow(float a, float b)
        {
            return b * b / (1 - a);
        }
        /// <summary>
        /// Implements the stamp function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Stamp(float a, float b)
        {
            return a + 2 * b - 1;
        }
        /// <summary>
        /// Implements the "freeze" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Freeze(float a, float b)
        {
            float x = 1 - a;
            return 1 - x * x / b;
        }
        /// <summary>
        /// Implements the "heat" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Heat(float a, float b)
        {
            float x = 1 - b;
            return 1 - x * x / a;
        }
        /// <summary>
        /// Implements the interpolation function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Interpolation(float a, float b)
        {
            return 0.5f - 0.25f * Maths.Cos(Maths.Pi * a) - 0.25f * Maths.Cos(Maths.Pi * b);
        }
        /// <summary>
        /// Implements the function of "soft light" (Adobe Photoshop).
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Photoshop(float a, float b)
        {
            if (b < 0.5)
            {
                return 2 * a * b + a * a * (1 - 2 * a);
            }
            return 2 * a * (1 - b) + Maths.Sqrt(a) * (2 * b - 1);
        }
        /// <summary>
        /// Implements the function of "soft light" (Illusions.hu).
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Illusions(float a, float b)
        {
            float x = 2 * (0.5f - b);
            float y = Maths.Pow(2, x);
            return Maths.Pow(a, y);
        }
        /// <summary>
        /// Implements the function of "soft light" (Pegtop).
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Pegtop(float a, float b)
        {
            return (1 - 2 * b) * a * a + 2 * b * a;
        }
        /// <summary>
        /// Implements the "Cairo" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Value</returns>
        public static float Fw3c(float a, float b)
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
        /// <param name="a">Value</param>
        /// <returns>Value</returns>
        public static float Gw3c(float a)
        {
            if (a <= 0.25)
            {
                return ((16 * a - 12) * a + 4) * a;
            }
            return Maths.Sqrt(a);
        }
        #endregion
    }
}
