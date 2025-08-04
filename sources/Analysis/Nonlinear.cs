using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements the solution of a nonlinear equation.
    /// <remarks>
    /// This class is a solution to the problem of finding the root of a nonlinear equation of the form F(x) = 0.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Nonlinear
    {
        #region Private data
        private Nonlinear.Method method;
        private float eps;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements the solution of a nonlinear equation.
        /// </summary>
        /// <param name="eps">Epsilon [0, 1]</param>
        /// <param name="method">Method for solving a nonlinear equation</param>
        public Nonlinear(float eps = 1e-8f, Nonlinear.Method method = Method.Secant)
        {
            this.method = method;
            this.Eps = eps;
        }
        /// <summary>
        /// Gets or sets the method for solving the nonlinear equation.
        /// </summary>
        public Nonlinear.Method MethodType
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
        /// Gets or sets the error value [0, 1].
        /// </summary>
        public float Eps
        {
            get
            {
                return this.eps;
            }
            set
            {
                this.eps = Maths.Float(value);
            }
        }
        /// <summary>
        /// Gets the root value of a nonlinear equation.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Start of line</param>
        /// <param name="b">End of line</param>
        /// <returns>Value</returns>
        public float Compute(IFloat function, float a, float b)
        {
            // chose method of nonlinear
            switch (method)
            {
                case Method.Chord:
                    return Nonlinear.chord(function, a, b, this.eps);
                case Method.FalsePosition:
                    return Nonlinear.falpo(function, a, b, this.eps);
                case Method.Secant:
                    return Nonlinear.secan(function, a, b, this.eps);

                default:
                    return Nonlinear.bisec(function, a, b, this.eps);
            }
        }
        /// <summary>
        /// Gets the root value of a nonlinear equation.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Start of line</param>
        /// <param name="b">End of line</param>
        /// <returns>Value</returns>
        public Complex32 Compute(IComplex function, Complex32 a, Complex32 b)
        {
            // chose method of nonlinear
            switch (method)
            {
                case Method.Chord:
                    return Nonlinear.chord(function, a, b, this.eps);
                case Method.FalsePosition:
                    return Nonlinear.falpo(function, a, b, this.eps);

                default:
                    return Nonlinear.secan(function, a, b, this.eps);
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
        /// <param name="eps"></param>
        /// <returns></returns>
        private static float bisec(IFloat f, float a, float b, float eps = 1e-8f)
        {
            float x1 = a; float x2 = b;
            float fb = f(b);
            float midpt;
            int n = 0;

            while (Math.Abs(x2 - x1) > eps && n < short.MaxValue)
            {
                midpt = 0.5f * (x1 + x2);

                if (fb * f(midpt) > 0)
                    x2 = midpt;
                else
                    x1 = midpt;
                n++;
            }
            return x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static float secan(IFloat f, float a, float b, float eps = 1e-8f)
        {
            float x1 = a;
            float x2 = b;
            float fb = f(b);
            float mpoint;
            int n = 0;

            while (Math.Abs(f(x2)) > eps && n < short.MaxValue)
            {
                mpoint = x2 - (x2 - x1) * fb / (fb - f(x1));
                x1 = x2;
                x2 = mpoint;
                fb = f(x2);
                n++;
            }
            return x2;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static float falpo(IFloat f, float a, float b, float eps = 1e-8f)
        {
            float x1 = a;
            float x2 = b;
            float fb = f(b);
            int n = 0;

            while (Math.Abs(x2 - x1) > eps && n < short.MaxValue)
            {
                float xpoint = x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
                if (fb * f(xpoint) > 0)
                    x2 = xpoint;
                else
                    x1 = xpoint;
                if (Math.Abs(f(xpoint)) < eps)
                    break;
                n++;
            }
            return x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static float chord(IFloat f, float a, float b, float eps = 1e-8f)
        {
            int n = 0;
            float x0 = (b - a) / 2.0f;
            float x;

            while (Math.Abs(f(x0) / b) > eps && n < short.MaxValue)
            {
                x = x0;
                x0 = x - (f(x) * (a - x)) / (f(a) - f(x));
                n++;
            }
            return x0;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static Complex32 chord(IComplex f, Complex32 a, Complex32 b, float eps = 1e-8f)
        {
            int n = 0;
            Complex32 x0 = (b - a) / 2.0;
            Complex32 x;

            while (Maths.Abs(f(x0) / b) > eps && n < short.MaxValue)
            {
                x = x0;
                x0 = x - (f(x) * (a - x)) / (f(a) - f(x));
                n++;
            }
            return x0;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static Complex32 secan(IComplex f, Complex32 a, Complex32 b, float eps = 1e-8f)
        {
            Complex32 x1 = a;
            Complex32 x2 = b;
            Complex32 fb = f(b);
            Complex32 mpoint;
            int n = 0;

            while (Maths.Abs(f(x2)) > eps && n < short.MaxValue)
            {
                mpoint = x2 - (x2 - x1) * fb / (fb - f(x1));
                x1 = x2;
                x2 = mpoint;
                fb = f(x2);
                n++;
            }
            return x2;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static Complex32 falpo(IComplex f, Complex32 a, Complex32 b, float eps = 1e-8f)
        {
            Complex32 x1 = a;
            Complex32 x2 = b;
            Complex32 fb = f(b);
            int n = 0;

            while (Maths.Abs(x2 - x1) > eps && n < short.MaxValue)
            {
                Complex32 xpoint = x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
                Complex32 fxpoint = f(xpoint);
                float s = fb.Real * fxpoint.Real;

                // sign
                if (s > 0)
                    x2 = xpoint;
                else
                    x1 = xpoint;

                if (Maths.Abs(fxpoint) < eps)
                    break;
                n++;
            }
            return x2 - (x2 - x1) * f(x2) / (f(x2) - f(x1));
        }
        #endregion

        #region Enums
        /// <summary>
        /// Method for solving a nonlinear equation.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Bisection method.
            /// </summary>
            Bisection,
            /// <summary>
            /// Chord method.
            /// </summary>
            Chord,
            /// <summary>
            /// Secant method.
            /// </summary>
            Secant,
            /// <summary>
            /// False position method.
            /// </summary>
            FalsePosition,
            #endregion
        }
        #endregion
    }
}
