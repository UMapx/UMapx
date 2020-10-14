using System;

namespace UMapx.Core
{
    /// <summary>
    /// Uses to work with kernel functions.
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
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Bicubic(double x)
        {
            if (x < 0)
            {
                x = -x;
            }

            double biCoef = 0;

            if (x <= 1)
            {
                biCoef = (1.5 * x - 2.5) * x * x + 1;
            }
            else if (x < 2)
            {
                biCoef = ((-0.5 * x + 2.5) * x - 4) * x + 2;
            }

            return biCoef;
        }
        #endregion

        #region Gaussian kernel function
        /// <summary>
        /// Returns the value of a Gaussian function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="sigma">Standard deviation (0, +inf)</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gaussian(double x, double sigma)
        {
            double t = x * x;
            double s = sigma * sigma;
            return Math.Exp(-t / s / 2.0);
        }
        /// <summary>
        /// Returns the value of a Gaussian function σ = 1.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gaussian(double x)
        {
            return Math.Exp(-x * x / 2);
        }
        #endregion

        #region Lanczos kernel function
        /// <summary>
        /// Returns the value of the Lanczos function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Parameter</param>
        /// <returns>Double precision floating point number</returns>
        public static double Lanczos(double x, double a)
        {
            if (x == 0)
            {
                return 1;
            }
            else if (-a <= x && x < a)
            {
                double pix = Math.PI * x;
                return a * Math.Sin(pix) * Math.Sin(pix / a) / (pix * pix);
            }
            return 0;
        }
        /// <summary>
        /// Returns the value of the Lanczos function, with a = 1.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Lanczos(double x)
        {
            return Lanczos(x, 1);
        }
        #endregion

        #region Uniform kernel function
        /// <summary>
        /// Returns the value of a uniform function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Uniform(double x)
        {
            double abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.5;
        }
        #endregion

        #region Triangular kernel function
        /// <summary>
        /// Returns the value of a triangular function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Triangular(double x)
        {
            double abs = Math.Abs(x);
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
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Trapezoid(double x)
        {
            double abs = Math.Abs(x);
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
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Epanechnikov(double x)
        {
            double abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.75 * (1 - x * x);
        }
        #endregion

        #region Quartic kernel function
        /// <summary>
        /// Returns the value of a Q function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Quartic(double x)
        {
            double abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.9375 * Math.Pow((1 - x * x), 2);
        }
        #endregion

        #region Triweight kernel function
        /// <summary>
        /// Returns the value of a T-function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Triweight(double x)
        {
            double abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 1.09375 * Math.Pow((1 - x * x), 3);
        }
        #endregion

        #region Tricube kernel function
        /// <summary>
        /// Returns the value of a tricubic function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Tricube(double x)
        {
            double abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.864197531 * Math.Pow((1 - x * x * x), 3);
        }
        #endregion

        #region Cosine kernel function
        /// <summary>
        /// Returns the value of the cosine function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Cosine(double x)
        {
            double abs = Math.Abs(x);
            if (abs > 1)
            {
                return 0;
            }
            return 0.7853981633975 * Math.Cos(1.570796326795 * x);
        }
        #endregion

        #region Logistic kernel function
        /// <summary>
        /// Returns the value of a logistic function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Logistic(double x)
        {
            return 1.0 / (Math.Exp(x) + 2 + Math.Exp(-x));
        }
        #endregion

        #region Sigmoid kernel function
        /// <summary>
        /// Returns the value of a sigmoid function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sigmoid(double x)
        {
            return 2.0 / Math.PI / (Math.Exp(x) + Math.Exp(-x));
        }
        #endregion

        #region Silverman function
        /// <summary>
        /// Returns the value of the Silverman function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Silverman(double x)
        {
            double abs = Math.Abs(x);
            double k = abs / 1.4142135623731;
            return 0.5 * Math.Exp(-k) * Math.Sin(k + 0.7853981633975);
        }
        #endregion
    }
}
