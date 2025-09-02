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
        private IntegrationMethod method;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements numerical integration.
        /// </summary>
        /// <param name="method">Integration method</param>
        public Integration(IntegrationMethod method = IntegrationMethod.Rectangle)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the integration method.
        /// </summary>
        public IntegrationMethod MethodType
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
                case IntegrationMethod.Midpoint:
                    return Integration.midp(function, a, b, n);

                case IntegrationMethod.Trapezoidal:
                    return Integration.trap(function, a, b, n);

                case IntegrationMethod.Simpson:
                    return Integration.simp(function, a, b, n);

                case IntegrationMethod.Romberg:
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
                case IntegrationMethod.Midpoint:
                    return Integration.midp(y, a, b, n);

                case IntegrationMethod.Trapezoidal:
                    return Integration.trap(y, a, b, n);

                case IntegrationMethod.Simpson:
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
        public Complex32 Compute(IComplex32 function, Complex32 a, Complex32 b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case IntegrationMethod.Midpoint:
                    return Integration.midp(function, a, b, n);

                case IntegrationMethod.Trapezoidal:
                    return Integration.trap(function, a, b, n);

                case IntegrationMethod.Simpson:
                    return Integration.simp(function, a, b, n);

                case IntegrationMethod.Romberg:
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
                case IntegrationMethod.Midpoint:
                    return Integration.midp(y, a, b, n);

                case IntegrationMethod.Trapezoidal:
                    return Integration.trap(y, a, b, n);

                case IntegrationMethod.Simpson:
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
            if (n < 2) return float.NaN;
            float h = (b - a) / n;
            float sum;

            if (n % 2 == 0)
            {
                float s_odd = 0f, s_even = 0f;
                for (int i = 1; i < n; i += 2) s_odd += f(a + i * h);
                for (int i = 2; i < n; i += 2) s_even += f(a + i * h);
                sum = (h / 3f) * (f(a) + f(b) + 4f * s_odd + 2f * s_even);
            }
            else
            {
                sum = (3f * h / 8f) * (f(a) + 3f * f(a + h) + 3f * f(a + 2 * h) + f(a + 3 * h));
                int m = n - 3;
                float s_odd = 0f, s_even = 0f;
                for (int j = 1; j < m; j += 2) s_odd += f(a + (3 + j) * h);
                for (int j = 2; j < m; j += 2) s_even += f(a + (3 + j) * h);
                sum += (h / 3f) * (f(a + 3 * h) + f(b) + 4f * s_odd + 2f * s_even);
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
            if (n < 3) return float.NaN;

            float h = (b - a) / (n - 1);
            float sum;

            if ((n - 1) % 2 == 0)
            {
                float s_odd = 0f, s_even = 0f;
                for (int i = 1; i < n - 1; i += 2) s_odd += y[i];
                for (int i = 2; i < n - 1; i += 2) s_even += y[i];
                sum = (h / 3f) * (y[0] + y[n - 1] + 4f * s_odd + 2f * s_even);
            }
            else
            {
                sum = (3f * h / 8f) * (y[0] + 3f * y[1] + 3f * y[2] + y[3]);
                float s_odd = 0f, s_even = 0f;
                for (int i = 4; i < n - 1; i += 2) s_odd += y[i - 1];
                for (int i = 5; i < n - 1; i += 2) s_even += y[i - 1];
                sum += (h / 3f) * (y[3] + y[n - 1] + 4f * s_odd + 2f * s_even);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="maxK"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static float romb(IFloat f, float a, float b, int maxK, float eps = 1e-8f)
        {
            if (maxK < 1) throw new ArgumentException();
            float[,] R = new float[maxK, maxK];

            float h = b - a;
            R[0, 0] = 0.5f * h * (f(a) + f(b));
            int n = 1;

            for (int k = 1; k < maxK; k++)
            {
                n *= 2; h *= 0.5f;
                float sum = 0f;
                for (int i = 1; i < n; i += 2) sum += f(a + i * h);
                R[k, 0] = 0.5f * R[k - 1, 0] + h * sum;

                float factor = 4f;
                for (int j = 1; j <= k; j++, factor *= 4f)
                    R[k, j] = (factor * R[k, j - 1] - R[k - 1, j - 1]) / (factor - 1f);

                if (Math.Abs(R[k, k] - R[k - 1, k - 1]) <= eps * Math.Abs(R[k, k]))
                    return R[k, k];
            }
            return R[maxK - 1, maxK - 1];
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        private static Complex32 rect(IComplex32 f, Complex32 a, Complex32 b, int n)
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
        private static Complex32 midp(IComplex32 f, Complex32 a, Complex32 b, int n)
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
        private static Complex32 trap(IComplex32 f, Complex32 a, Complex32 b, int n)
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
        private static Complex32 simp(IComplex32 f, Complex32 a, Complex32 b, int n)
        {
            if (n < 2) return Complex32.NaN;
            Complex32 h = (b - a) / n;
            Complex32 sum;

            if (n % 2 == 0)
            {
                Complex32 s_odd = 0f, s_even = 0f;
                for (int i = 1; i < n; i += 2) s_odd += f(a + i * h);
                for (int i = 2; i < n; i += 2) s_even += f(a + i * h);
                sum = (h / 3f) * (f(a) + f(b) + 4f * s_odd + 2f * s_even);
            }
            else
            {
                sum = (3f * h / 8f) * (f(a) + 3f * f(a + h) + 3f * f(a + 2 * h) + f(a + 3 * h));
                int m = n - 3;
                Complex32 s_odd = 0f, s_even = 0f;
                for (int j = 1; j < m; j += 2) s_odd += f(a + (3 + j) * h);
                for (int j = 2; j < m; j += 2) s_even += f(a + (3 + j) * h);
                sum += (h / 3f) * (f(a + 3 * h) + f(b) + 4f * s_odd + 2f * s_even);
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
            if (n < 3) return Complex32.NaN;

            Complex32 h = (b - a) / (n - 1);
            Complex32 sum;

            if ((n - 1) % 2 == 0)
            {
                Complex32 s_odd = 0f, s_even = 0f;
                for (int i = 1; i < n - 1; i += 2) s_odd += y[i];
                for (int i = 2; i < n - 1; i += 2) s_even += y[i];
                sum = (h / 3f) * (y[0] + y[n - 1] + 4f * s_odd + 2f * s_even);
            }
            else
            {
                sum = (3f * h / 8f) * (y[0] + 3f * y[1] + 3f * y[2] + y[3]);
                Complex32 s_odd = 0f, s_even = 0f;
                for (int i = 4; i < n - 1; i += 2) s_odd += y[i - 1];
                for (int i = 5; i < n - 1; i += 2) s_even += y[i - 1];
                sum += (h / 3f) * (y[3] + y[n - 1] + 4f * s_odd + 2f * s_even);
            }
            return sum;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="maxK"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static Complex32 romb(IComplex32 f, Complex32 a, Complex32 b, int maxK, float eps = 1e-8f)
        {
            if (maxK < 1) throw new ArgumentException();
            Complex32[,] R = new Complex32[maxK, maxK];

            Complex32 h = b - a;
            R[0, 0] = 0.5f * h * (f(a) + f(b));
            int n = 1;

            for (int k = 1; k < maxK; k++)
            {
                n *= 2; h *= 0.5f;
                Complex32 sum = 0f;
                for (int i = 1; i < n; i += 2) sum += f(a + i * h);
                R[k, 0] = 0.5f * R[k - 1, 0] + h * sum;

                Complex32 factor = 4f;
                for (int j = 1; j <= k; j++, factor *= 4f)
                    R[k, j] = (factor * R[k, j - 1] - R[k - 1, j - 1]) / (factor - 1f);

                if (Maths.Abs(R[k, k] - R[k - 1, k - 1]) <= eps * Maths.Abs(R[k, k]))
                    return R[k, k];
            }
            return R[maxK - 1, maxK - 1];
        }
        #endregion
    }
}
