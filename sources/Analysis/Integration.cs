using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements numerical integration.
    /// <remarks>
    /// This class is a solution to the problem of finding the value of the integral of the function F(x) within the values of a and b.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Integration
    {
        #region Private data
        private Integration.Method method;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements numerical integration.
        /// </summary>
        /// <param name="method">Integration method</param>
        public Integration(Integration.Method method = Method.Rectangle)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the integration method.
        /// </summary>
        public Integration.Method MethodType
        {
            get
            {
                return this.method;
            }
            set
            {
                this.method = value;
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Value</returns>
        public float Compute(IFloat function, float a, float b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case Method.Midpoint:
                    return Integration.midp(function, a, b, n);

                case Method.Trapezoidal:
                    return Integration.trap(function, a, b, n);

                case Method.Simpson:
                    return Integration.simp(function, a, b, n);

                case Method.Romberg:
                    return Integration.romb(function, a, b, n);

                default:
                    return Integration.rect(function, a, b, n);
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="y">Function vector</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Value</returns>
        public float Compute(float[] y, float a, float b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case Method.Midpoint:
                    return Integration.midp(y, a, b, n);

                case Method.Trapezoidal:
                    return Integration.trap(y, a, b, n);

                case Method.Simpson:
                    return Integration.simp(y, a, b, n);

                default:
                    return Integration.rect(y, a, b, n);
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Complex number</returns>
        public Complex32 Compute(IComplex function, Complex32 a, Complex32 b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case Method.Midpoint:
                    return Integration.midp(function, a, b, n);

                case Method.Trapezoidal:
                    return Integration.trap(function, a, b, n);

                case Method.Simpson:
                    return Integration.simp(function, a, b, n);

                case Method.Romberg:
                    return Integration.romb(function, a, b, n);

                default:
                    return Integration.rect(function, a, b, n);
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="y">Function vector</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Complex number</returns>
        public Complex32 Compute(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case Method.Midpoint:
                    return Integration.midp(y, a, b, n);

                case Method.Trapezoidal:
                    return Integration.trap(y, a, b, n);

                case Method.Simpson:
                    return Integration.simp(y, a, b, n);

                default:
                    return Integration.rect(y, a, b, n);
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static float rect(IFloat f, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + i * h);
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static float rect(float[] y, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * y[i];
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static float midp(IFloat f, float a, float b, int n)
        {
            // Midpoint
            float sum = 0.0f;
            float h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + (i + 0.5f) * h);
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static float midp(float[] y, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += h * 0.5f * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static float trap(IFloat f, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += 0.5f * h * (f(a + i * h) + f(a + (i + 1) * h));
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static float trap(float[] y, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += 0.5f * h * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static float simp(IFloat f, float a, float b, int n)
        {
            if (n < 3) return float.NaN; //Need at least 3 points
            float sum = 0.0f;
            float h = (b - a) / n;
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 1; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            else
            {
                sum = 3 * h * (f(a) + 3 * f(a + h) + 3 * f(a + 2 * h) + f(a + 3 * h)) / 8;
                for (int i = 3; i < n - 1; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static float simp(float[] y, float a, float b, int n)
        {
            float h = (b - a) / n;
            //Need at least 3 points
            if (n < 3 || h == 0) return float.NaN;
            float sum = 0.0f;
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 1; i += 2)
                {
                    sum += h * (y[i] + 4 * y[i + 1] + y[i + 2]) / 3;
                }
            }
            else
            {
                sum = 3 * h * (y[0] + 3 * y[1] + 3 * y[2] + y[3]) / 8;
                for (int i = 3; i < n - 1; i += 2)
                {
                    sum += h * (y[i] + 4 * y[i + 1] + y[i + 2]) / 3;
                }
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="iterations"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static float romb(IFloat f, float a, float b, int iterations, float eps = 1e-8f)
        {
            int n = 2;
            float h = b - a;
            float sum = 0.0f;
            int j = 0;
            float[,] R = new float[iterations, iterations];
            R[1, 1] = h * (f(a) + f(b)) / 2.0f;
            h = h / 2;
            R[2, 1] = R[1, 1] / 2 + h * f(a + h);
            R[2, 2] = (4 * R[2, 1] - R[1, 1]) / 3;
            for (j = 3; j <= iterations; j++)
            {
                n = 2 * n;
                h = h / 2;
                sum = 0.0f;
                for (int k = 1; k <= n; k += 2)
                {
                    sum += f(a + k * h);
                }
                R[j, 1] = R[j - 1, 1] / 2 + h * sum;
                float factor = 4.0f;
                for (int k = 2; k <= j; k++)
                {
                    R[j, k] = (factor * R[j, k - 1] - R[j - 1, k - 1]) / (factor - 1);
                    factor = factor * 4.0f;
                }
                if (Math.Abs(R[j, j] - R[j, j - 1]) < eps * Math.Abs(R[j, j]))
                {
                    sum = R[j, j];
                    return sum;
                }
            }
            sum = R[n, n];
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex32 rect(IComplex f, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + i * h);
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex32 rect(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * y[i];
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex32 midp(IComplex f, Complex32 a, Complex32 b, int n)
        {
            // Midpoint
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + (i + 0.5) * h);
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex32 midp(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += h * 0.5 * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex32 trap(IComplex f, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += 0.5 * h * (f(a + i * h) + f(a + (i + 1) * h));
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex32 trap(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += 0.5 * h * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex32 simp(IComplex f, Complex32 a, Complex32 b, int n)
        {
            if (n < 3) return float.NaN; //Need at least 3 points
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / n;
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 1; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            else
            {
                sum = 3 * h * (f(a) + 3 * f(a + h) + 3 * f(a + 2 * h) + f(a + 3 * h)) / 8;
                for (int i = 3; i < n - 1; i += 2)
                {
                    sum += h * (f(a + i * h) + 4 * f(a + (i + 1) * h) + f(a + (i + 2) * h)) / 3;
                }
            }
            return sum;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex32 simp(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            Complex32 h = (b - a) / n;
            //Need at least 3 points
            if (n < 3 || h == 0) return float.NaN;
            Complex32 sum = 0.0;
            if (n % 2 != 0)
            {
                for (int i = 0; i < n - 1; i += 2)
                {
                    sum += h * (y[i] + 4 * y[i + 1] + y[i + 2]) / 3;
                }
            }
            else
            {
                sum = 3 * h * (y[0] + 3 * y[1] + 3 * y[2] + y[3]) / 8;
                for (int i = 3; i < n - 1; i += 2)
                {
                    sum += h * (y[i] + 4 * y[i + 1] + y[i + 2]) / 3;
                }
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="iterations"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static Complex32 romb(IComplex f, Complex32 a, Complex32 b, int iterations, float eps = 1e-8f)
        {
            int n = 2;
            Complex32 h = b - a;
            Complex32 sum = 0.0;
            int j = 0;
            Complex32[,] R = new Complex32[iterations, iterations];
            R[1, 1] = h * (f(a) + f(b)) / 2.0;
            h = h / 2;
            R[2, 1] = R[1, 1] / 2 + h * f(a + h);
            R[2, 2] = (4 * R[2, 1] - R[1, 1]) / 3;
            for (j = 3; j <= iterations; j++)
            {
                n = 2 * n;
                h = h / 2;
                sum = 0.0;
                for (int k = 1; k <= n; k += 2)
                {
                    sum += f(a + k * h);
                }
                R[j, 1] = R[j - 1, 1] / 2 + h * sum;
                float factor = 4.0f;
                for (int k = 2; k <= j; k++)
                {
                    R[j, k] = (factor * R[j, k - 1] - R[j - 1, k - 1]) / (factor - 1);
                    factor = factor * 4.0f;
                }
                if (Maths.Abs(R[j, j] - R[j, j - 1]) < eps * Maths.Abs(R[j, j]))
                {
                    sum = R[j, j];
                    return sum;
                }
            }
            sum = R[n, n];
            return sum;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Integration method.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Rectangle method.
            /// </summary>
            Rectangle,
            /// <summary>
            /// Midpoint method.
            /// </summary>
            Midpoint,
            /// <summary>
            /// Trapezoidal method.
            /// </summary>
            Trapezoidal,
            /// <summary>
            /// Simpson method.
            /// </summary>
            Simpson,
            /// <summary>
            /// Romberg method.
            /// </summary>
            Romberg,
            #endregion
        }
        #endregion
    }
}
