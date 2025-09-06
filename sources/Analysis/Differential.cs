using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements a solution to a differential equation.
    /// </summary>
    /// <remarks>
    /// This class is a solution to the Cauchy problem for the ordinary differential equation y' = F(x, y).
    /// </remarks>
    [Serializable]
    public class Differential
    {
        #region Private data
        private DifferentialMethod method;
        #endregion

        #region Differentiation components
        /// <summary>
        /// Initializes a class that implements the solution of a differential equation.
        /// </summary>
        /// <param name="method">Differentiation method</param>
        public Differential(DifferentialMethod method = DifferentialMethod.RungeKutta4)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the differentiation method.
        /// </summary>
        public DifferentialMethod MethodType
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
            // choose method of differentiation
            switch (method)
            {
                case DifferentialMethod.Euler:
                    return Differential.Euler(function, x, y0);

                case DifferentialMethod.Fehlberg:
                    return Differential.Fehlberg(function, x, y0);

                case DifferentialMethod.RungeKutta4:
                    return Differential.RungeKutta4(function, x, y0);

                default:
                    return Differential.RungeKutta2(function, x, y0);
            }
        }
        /// <summary>
        /// Returns the value of a differential equation.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <returns>Array of function values</returns>
        public Complex32[] Compute(IComplex32Mesh function, Complex32[] x, Complex32 y0)
        {
            // choose method of differentiation
            switch (method)
            {
                case DifferentialMethod.Euler:
                    return Differential.Euler(function, x, y0);

                case DifferentialMethod.Fehlberg:
                    return Differential.Fehlberg(function, x, y0);

                case DifferentialMethod.RungeKutta4:
                    return Differential.RungeKutta4(function, x, y0);

                default:
                    return Differential.RungeKutta2(function, x, y0);
            }
        }
        #endregion

        #region Recompute voids
        /// <summary>
        /// Returns the value of a differential equation calculated by the Adams-Bashforth method.
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
                    float hnext = x[i + 1] - x[i];
                    float sum = y[i - 1];

                    for (j = 0; j < order; j++)
                    {
                        float t = x[i - j];
                        sum += hnext * c[j] * function(t, y[i - j - 1]);
                    }

                    y[i] = sum;
                }

                return y;
            }

            // classic differential
            return this.Compute(function, x, y0);
        }
        /// <summary>
        /// Returns the value of a differential equation calculated by the Adams-Bashforth method.
        /// </summary>
        /// <param name="function">The delegate of a continuous function depending on two variables</param>
        /// <param name="x">Array of values argument</param>
        /// <param name="y0">Value</param>
        /// <param name="order">Order</param>
        /// <returns>Array of function values</returns>
        public Complex32[] Compute(IComplex32Mesh function, Complex32[] x, Complex32 y0, int order = 2)
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
                    Complex32 hnext = x[i + 1] - x[i];
                    Complex32 sum = y[i - 1];

                    for (j = 0; j < order; j++)
                    {
                        Complex32 t = x[i - j];
                        sum += hnext * c[j] * function(t, y[i - j - 1]);
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
        /// Returns an array of coefficient values for the Adams-Bashforth formula.
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
                    A[i, j] = Maths.Pow(j, i);
                }
                c[i] = Maths.Pow(-1, i) / (i + 1);
            }

            return A.Solve(c);
        }
        #endregion

        #region Runge-Kutta
        /// <summary>
        /// Explicit Euler method (1st order) for y' = f(x, y) on a given grid x.
        /// </summary>
        /// <remarks>
        /// Works with nonuniform grids. The returned array has length n = x.Length - 1
        /// and contains y at x[1], x[2], ..., x[n]. The initial value y0 (at x[0]) is not included.
        /// </remarks>
        /// <param name="f">Right-hand side f(x, y)</param>
        /// <param name="x">Monotone grid points</param>
        /// <param name="y0">Initial value y(x[0])</param>
        /// <returns>Solution values at x[1..n]</returns>
        private static float[] Euler(IFloatMesh f, float[] x, float y0)
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
        /// Runge–Kutta method of order 2 (midpoint form) on grid x.
        /// </summary>
        /// <remarks>
        /// Works with nonuniform grids. Returns y at x[1..n]; y0 is not included.
        /// </remarks>
        /// <param name="f">Right-hand side f(x, y)</param>
        /// <param name="x">Grid points</param>
        /// <param name="y0">Initial value y(x[0])</param>
        /// <returns>Solution values at x[1..n]</returns>
        private static float[] RungeKutta2(IFloatMesh f, float[] x, float y0)
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
        /// Classic Runge–Kutta method of order 4 (RK4) on grid x.
        /// </summary>
        /// <remarks>
        /// Works with nonuniform grids. Returns y at x[1..n]; y0 is not included.
        /// </remarks>
        /// <param name="f">Right-hand side f(x, y)</param>
        /// <param name="x">Grid points</param>
        /// <param name="y0">Initial value y(x[0])</param>
        /// <returns>Solution values at x[1..n]</returns>
        private static float[] RungeKutta4(IFloatMesh f, float[] x, float y0)
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
        /// Fehlberg (RKF45 scheme without stepsize adaptation): fixed-step 4th-order estimate.
        /// </summary>
        /// <remarks>
        /// Uses the embedded Runge–Kutta-Fehlberg coefficients on a fixed step.
        /// No local error return or adaptive control. Returns y at x[1..n]; y0 is not included.
        /// </remarks>
        /// <param name="f">Right-hand side f(x, y)</param>
        /// <param name="x">Grid points</param>
        /// <param name="y0">Initial value y(x[0])</param>
        /// <returns>Solution values at x[1..n]</returns>
        private static float[] Fehlberg(IFloatMesh f, float[] x, float y0)
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
        /// Explicit Euler method (1st order) for y' = f(x, y) on a given grid x.
        /// </summary>
        /// <remarks>
        /// Works with nonuniform grids. The returned array has length n = x.Length - 1
        /// and contains y at x[1], x[2], ..., x[n]. The initial value y0 (at x[0]) is not included.
        /// </remarks>
        /// <param name="f">Right-hand side f(x, y)</param>
        /// <param name="x">Monotone grid points</param>
        /// <param name="y0">Initial value y(x[0])</param>
        /// <returns>Solution values at x[1..n]</returns>
        private static Complex32[] Euler(IComplex32Mesh f, Complex32[] x, Complex32 y0)
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
        /// Runge–Kutta method of order 2 (midpoint form) on grid x.
        /// </summary>
        /// <remarks>
        /// Works with nonuniform grids. Returns y at x[1..n]; y0 is not included.
        /// </remarks>
        /// <param name="f">Right-hand side f(x, y)</param>
        /// <param name="x">Grid points</param>
        /// <param name="y0">Initial value y(x[0])</param>
        /// <returns>Solution values at x[1..n]</returns>
        private static Complex32[] RungeKutta2(IComplex32Mesh f, Complex32[] x, Complex32 y0)
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
        /// Classic Runge–Kutta method of order 4 (RK4) on grid x.
        /// </summary>
        /// <remarks>
        /// Works with nonuniform grids. Returns y at x[1..n]; y0 is not included.
        /// </remarks>
        /// <param name="f">Right-hand side f(x, y)</param>
        /// <param name="x">Grid points</param>
        /// <param name="y0">Initial value y(x[0])</param>
        /// <returns>Solution values at x[1..n]</returns>
        private static Complex32[] RungeKutta4(IComplex32Mesh f, Complex32[] x, Complex32 y0)
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
        /// Fehlberg (RKF45 scheme without stepsize adaptation): fixed-step 4th-order estimate.
        /// </summary>
        /// <remarks>
        /// Uses the embedded Runge–Kutta-Fehlberg coefficients on a fixed step.
        /// No local error return or adaptive control. Returns y at x[1..n]; y0 is not included.
        /// </remarks>
        /// <param name="f">Right-hand side f(x, y)</param>
        /// <param name="x">Grid points</param>
        /// <param name="y0">Initial value y(x[0])</param>
        /// <returns>Solution values at x[1..n]</returns>
        private static Complex32[] Fehlberg(IComplex32Mesh f, Complex32[] x, Complex32 y0)
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
    }
}
