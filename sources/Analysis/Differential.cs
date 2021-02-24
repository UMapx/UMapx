using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements a solution to a differential equation.
    /// <remarks>
    /// This class is a solution to the Cauchy problem for the ordinary differential equation y' = F(x, y).
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Differential
    {
        #region Private data
        private Differential.Method method;
        #endregion

        #region Diferentiation components
        /// <summary>
        /// Initializes a class that implements the solution of a differential equation.
        /// </summary>
        /// <param name="method">Differentiation method</param>
        public Differential(Differential.Method method = Method.RungeKutta4)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the differentiation method.
        /// </summary>
        public Differential.Method MethodType
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
        /// Returns the value of a differential equation.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <returns>Array of function values</returns>
        public float[] Compute(IFloatMesh function, float[] x, float y0)
        {
            // chose method of differentiation
            switch (method)
            {
                case Method.Euler:
                    return Differential.euler(function, x, y0);

                case Method.Fehlberg:
                    return Differential.fehlberg(function, x, y0);

                case Method.RungeKutta4:
                    return Differential.rungeKutta4(function, x, y0);

                default:
                    return Differential.rungeKutta2(function, x, y0);
            }
        }
        /// <summary>
        /// Returns the value of a differential equation.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <returns>Array of function values</returns>
        public Complex32[] Compute(IComplexMesh function, Complex32[] x, Complex32 y0)
        {
            // chose method of differentiation
            switch (method)
            {
                case Method.Euler:
                    return Differential.euler(function, x, y0);

                case Method.Fehlberg:
                    return Differential.fehlberg(function, x, y0);

                case Method.RungeKutta4:
                    return Differential.rungeKutta4(function, x, y0);

                default:
                    return Differential.rungeKutta2(function, x, y0);
            }
        }
        #endregion

        #region Recompute voids
        /// <summary>
        /// Returns the value of a differential equation calculated by the Adams-Bashfort method.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <param name="order">Order</param>
        /// <returns>Array of function values</returns>
        public float[] Compute(IFloatMesh function, float[] x, float y0, int order = 2)
        {
            int n = x.Length - 1;

            // if order more than 1
            // Adams-Bashfort method
            if (order > 1 && order < n)
            {
                // params
                int i, j, k = order + 1;
                float[] y = new float[n];
                float[] r = new float[k];
                float[] c = Differential.GetCoefficients(order);
                float h, t, sum;

                // compute first points by order
                for (i = 0; i < k; i++)
                    r[i] = x[i];

                // classic differential
                r = this.Compute(function, r, y0);

                for (i = 0; i < order; i++)
                    y[i] = r[i];

                // Adams-Bashforth method
                // for order
                for (i = order; i < n; i++)
                {
                    sum = y[i - 1];

                    for (j = 0; j < order; j++)
                    {
                        t = x[i - j];
                        h = t - x[i - j - 1];
                        sum += h * c[j] * function(t, y[i - j - 1]);
                    }

                    y[i] = sum;
                }

                return y;
            }

            // classic differential
            return this.Compute(function, x, y0);
        }
        /// <summary>
        /// Returns the value of a differential equation calculated by the Adams-Bashfort method.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <param name="order">Order</param>
        /// <returns>Array of function values</returns>
        public Complex32[] Compute(IComplexMesh function, Complex32[] x, Complex32 y0, int order = 2)
        {
            int n = x.Length - 1;

            // if order more than 1
            // Adams-Bashfort method
            if (order > 1 && order < n)
            {
                // params
                int i, j, k = order + 1;
                Complex32[] y = new Complex32[n];
                Complex32[] r = new Complex32[k];
                float[] c = Differential.GetCoefficients(order);
                Complex32 h, t, sum;

                // compute first points by order
                for (i = 0; i < k; i++)
                    r[i] = x[i];

                // classic differential
                r = this.Compute(function, r, y0);

                for (i = 0; i < order; i++)
                    y[i] = r[i];

                // Adams-Bashforth method
                // for order
                for (i = order; i < n; i++)
                {
                    sum = y[i - 1];

                    for (j = 0; j < order; j++)
                    {
                        t = x[i - j];
                        h = t - x[i - j - 1];
                        sum += h * c[j] * function(t, y[i - j - 1]);
                    }

                    y[i] = sum;
                }

                return y;
            }

            // classic differential
            return this.Compute(function, x, y0);
        }
        #endregion

        #region Adams-Bashforth
        /// <summary>
        /// Returns an array of coefficient values for the Adams-Bashfort formula.
        /// </summary>
        /// <param name="order">Order</param>
        /// <returns>Array</returns>
        public static float[] GetCoefficients(int order)
        {
            float[,] A = new float[order, order];
            float[] c = new float[order];
            int i, j;

            for (i = 0; i < order; i++)
            {
                for (j = 0; j < order; j++)
                {
                    A[i, j] = (float)Math.Pow(j, i);
                }
                c[i] = (float)Math.Pow(-1, i) / (i + 1);
            }

            return A.Solve(c);
        }
        #endregion

        #region Runge-Kutta
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static float[] euler(IFloatMesh f, float[] x, float y0)
        {
            int n = x.Length - 1;
            float xnew, ynew = y0, h;
            float[] result = new float[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                ynew = ynew + f(xnew, ynew) * h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static float[] rungeKutta2(IFloatMesh f, float[] x, float y0)
        {
            int n = x.Length - 1;
            float xnew, ynew = y0, h, k1, k2;
            float[] result = new float[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.5f * h, ynew + 0.5f * k1);
                ynew = ynew + k2;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static float[] rungeKutta4(IFloatMesh f, float[] x, float y0)
        {
            int n = x.Length - 1;
            float xnew, ynew = y0, h, k1, k2, k3, k4;
            float[] result = new float[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.5f * h, ynew + 0.5f * k1);
                k3 = h * f(xnew + 0.5f * h, ynew + 0.5f * k2);
                k4 = h * f(xnew + h, ynew + k3);
                ynew = ynew + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static float[] fehlberg(IFloatMesh f, float[] x, float y0)
        {
            int n = x.Length - 1;
            float xnew, ynew = y0, h, k1, k2, k3, k4, k5, k6;
            float[] result = new float[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.25f * h, ynew + 0.25f * k1);
                k3 = h * f(xnew + 3 * h / 8, ynew + 3 * k1 / 32 + 9 * k2 / 32);
                k4 = h * f(xnew + 12 * h / 13, ynew + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197);
                k5 = h * f(xnew + h, ynew + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104);
                k6 = h * f(xnew + 0.5f * h, ynew - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40);
                ynew = ynew + 25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - 0.2f * k5;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static Complex32[] euler(IComplexMesh f, Complex32[] x, Complex32 y0)
        {
            int n = x.Length - 1;
            Complex32 xnew, ynew = y0, h;
            Complex32[] result = new Complex32[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                ynew = ynew + f(xnew, ynew) * h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static Complex32[] rungeKutta2(IComplexMesh f, Complex32[] x, Complex32 y0)
        {
            int n = x.Length - 1;
            Complex32 xnew, ynew = y0, h, k1, k2;
            Complex32[] result = new Complex32[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.5 * h, ynew + 0.5 * k1);
                ynew = ynew + k2;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static Complex32[] rungeKutta4(IComplexMesh f, Complex32[] x, Complex32 y0)
        {
            int n = x.Length - 1;
            Complex32 xnew, ynew = y0, h, k1, k2, k3, k4;
            Complex32[] result = new Complex32[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.5 * h, ynew + 0.5 * k1);
                k3 = h * f(xnew + 0.5 * h, ynew + 0.5 * k2);
                k4 = h * f(xnew + h, ynew + k3);
                ynew = ynew + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="x"></param>
        /// <param name="y0"></param>
        /// <returns></returns>
        private static Complex32[] fehlberg(IComplexMesh f, Complex32[] x, Complex32 y0)
        {
            int n = x.Length - 1;
            Complex32 xnew, ynew = y0, h, k1, k2, k3, k4, k5, k6;
            Complex32[] result = new Complex32[n];

            for (int i = 0; i < n; i++)
            {
                h = x[i + 1] - x[i];
                xnew = x[i];
                k1 = h * f(xnew, ynew);
                k2 = h * f(xnew + 0.25 * h, ynew + 0.25 * k1);
                k3 = h * f(xnew + 3 * h / 8, ynew + 3 * k1 / 32 + 9 * k2 / 32);
                k4 = h * f(xnew + 12 * h / 13, ynew + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197);
                k5 = h * f(xnew + h, ynew + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104);
                k6 = h * f(xnew + 0.5 * h, ynew - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40);
                ynew = ynew + 25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - 0.2 * k5;
                xnew = xnew + h;
                result[i] = ynew;
            }
            return result;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Differentiation method
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Euler method.
            /// </summary>
            Euler,
            /// <summary>
            /// The second-order Runge-Kutta method.
            /// </summary>
            RungeKutta2,
            /// <summary>
            /// Fourth-order Runge-Kutta method.
            /// </summary>
            RungeKutta4,
            /// <summary>
            /// Felberg's method.
            /// </summary>
            Fehlberg,
            #endregion
        }
        #endregion
    }
}
