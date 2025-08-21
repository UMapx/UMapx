using System;

namespace UMapx.Core
{
    /// <summary>
    /// Used to work with kernel functions.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Kernel_(statistics)
    /// </remarks>
    /// </summary>
    public static class Kernel
    {
        #region Bicubic kernel function
        /// <summary>
        /// Returns the value of a bicubic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Bicubic(float x)
        {
            if (x < 0)
            {
                x = -x;
            }

            float biCoef = 0;

            if (x <= 1)
            {
                biCoef = (1.5f * x - 2.5f) * x * x + 1;
            }
            else if (x < 2)
            {
                biCoef = ((-0.5f * x + 2.5f) * x - 4) * x + 2;
            }

            return biCoef;
        }
        #endregion

        #region Gaussian kernel function
        /// <summary>
        /// Returns the value of a Gaussian function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="sigma">Standard deviation (0, +inf)</param>
        /// <returns>Value</returns>
        public static float Gaussian(float x, float sigma)
        {
            float t = x * x;
            float s = sigma * sigma;
            return (float)Math.Exp(-t / s / 2.0);
        }
        /// <summary>
        /// Returns the value of a Gaussian function σ = 1.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Gaussian(float x)
        {
            return (float)Math.Exp(-x * x / 2);
        }
        #endregion

        #region Lanczos kernel function
        /// <summary>
        /// Returns the value of the Lanczos function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="a">Parameter</param>
        /// <returns>Value</returns>
        public static float Lanczos(float x, float a)
        {
            if (x == 0)
            {
                return 1;
            }
            else if (-a <= x && x < a)
            {
                float pix = Maths.Pi * x;
                return a * (float)Math.Sin(pix) * (float)Math.Sin(pix / a) / (pix * pix);
            }
            return 0;
        }
        /// <summary>
        /// Returns the value of the Lanczos function, with a = 1.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Lanczos(float x)
        {
            return Lanczos(x, 1);
        }
        #endregion

        #region Uniform kernel function
        /// <summary>
        /// Returns the value of a uniform function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Uniform(float x)
        {
            float abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.5f;
        }
        #endregion

        #region Triangular kernel function
        /// <summary>
        /// Returns the value of a triangular function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Triangular(float x)
        {
            float abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 1 - abs;
        }
        #endregion

        #region Trapezoid kernel function
        /// <summary>
        /// Returns the value of the trapezoid function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Trapezoid(float x)
        {
            float abs = Math.Abs(x);
            if (abs < 1.0 / 2)
            {
                return 1;
            }
            else if (abs < 1.0)
            {
                return 2 * (1 - abs);
            }
            return 0;
        }
        #endregion

        #region Epanechnikov kernel function
        /// <summary>
        /// Returns the value of the Epanechnikov function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Epanechnikov(float x)
        {
            float abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.75f * (1 - x * x);
        }
        #endregion

        #region Quartic kernel function
        /// <summary>
        /// Returns the value of a Q function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Quartic(float x)
        {
            float abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.9375f * (float)Math.Pow((1 - x * x), 2);
        }
        #endregion

        #region Triweight kernel function
        /// <summary>
        /// Returns the value of a T-function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Triweight(float x)
        {
            float abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 1.09375f * (float)Math.Pow((1 - x * x), 3);
        }
        #endregion

        #region Tricube kernel function
        /// <summary>
        /// Returns the value of a tricubic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Tricube(float x)
        {
            float abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.864197531f * (float)Math.Pow((1 - x * x * x), 3);
        }
        #endregion

        #region Cosine kernel function
        /// <summary>
        /// Returns the value of the cosine function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Cosine(float x)
        {
            float abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.7853981633975f * (float)Math.Cos(1.570796326795 * x);
        }
        #endregion

        #region Logistic kernel function
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Logistic(float x)
        {
            return 1.0f / (float)(Math.Exp(x) + 2 + Math.Exp(-x));
        }
        #endregion

        #region Sigmoid kernel function
        /// <summary>
        /// Returns the value of a sigmoid function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Sigmoid(float x)
        {
            return 2.0f / Maths.Pi / (float)(Math.Exp(x) + Math.Exp(-x));
        }
        #endregion

        #region Silverman function
        /// <summary>
        /// Returns the value of the Silverman function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public static float Silverman(float x)
        {
            float abs = Math.Abs(x);
            float k = abs / 1.4142135623731f;
            return 0.5f * (float)Math.Exp(-k) * (float)Math.Sin(k + 0.7853981633975);
        }
        #endregion
    }
}
