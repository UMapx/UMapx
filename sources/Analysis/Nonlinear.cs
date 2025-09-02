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
        private NonlinearMethod method;
        private float eps;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements the solution of a nonlinear equation.
        /// </summary>
        /// <param name="eps">Epsilon [0, 1]</param>
        /// <param name="method">Method for solving a nonlinear equation</param>
        public Nonlinear(float eps = 1e-8f, NonlinearMethod method = NonlinearMethod.Secant)
        {
            this.method = method;
            this.Eps = eps;
        }
        /// <summary>
        /// Gets or sets the method for solving the nonlinear equation.
        /// </summary>
        public NonlinearMethod MethodType
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
                case NonlinearMethod.Chord:
                    return Nonlinear.chord(function, a, b, this.eps);
                case NonlinearMethod.FalsePosition:
                    return Nonlinear.falpo(function, a, b, this.eps);
                case NonlinearMethod.Secant:
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
        public Complex32 Compute(IComplex32 function, Complex32 a, Complex32 b)
        {
            // chose method of nonlinear
            switch (method)
            {
                case NonlinearMethod.Chord:
                    return Nonlinear.chord(function, a, b, this.eps);
                case NonlinearMethod.FalsePosition:
                    throw new NotSupportedException("False position is not defined for complex-valued functions");

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
            float x1 = a, x2 = b;
            float fa = f(x1), fb = f(x2);
            if (fa == 0f) return x1;
            if (fb == 0f) return x2;
            if (fa * fb > 0f) throw new ArgumentException("Bisection requires f(a)*f(b) <= 0");

            int n = 0;
            while ((Math.Abs(x2 - x1) > eps) && (n++ < short.MaxValue))
            {
                float mid = 0.5f * (x1 + x2);
                float fm = f(mid);
                if (fm == 0f) return mid;

                if (fa * fm < 0f) { x2 = mid; fb = fm; }
                else { x1 = mid; fa = fm; }
            }
            float denom = (f(x2) - f(x1));
            return denom != 0f ? x2 - (x2 - x1) * f(x2) / denom
                               : 0.5f * (x1 + x2);
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
            float x1 = a, x2 = b;
            float f1 = f(x1), f2 = f(x2);
            int n = 0;

            while (Math.Abs(f2) > eps && Math.Abs(x2 - x1) > eps && n++ < short.MaxValue)
            {
                float denom = (f2 - f1);
                if (denom == 0f) break;
                float x3 = x2 - (x2 - x1) * f2 / denom;

                x1 = x2; f1 = f2;
                x2 = x3; f2 = f(x2);
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
            float x1 = a, x2 = b;
            float fa = f(x1), fb = f(x2);
            if (fa == 0f) return x1;
            if (fb == 0f) return x2;
            if (fa * fb > 0f) throw new ArgumentException("False position requires f(a)*f(b) <= 0");

            int n = 0;
            float x = x2;
            while (n++ < short.MaxValue)
            {
                float denom = (fb - fa);
                if (denom == 0f) break;
                x = x2 - (x2 - x1) * fb / denom;
                float fx = f(x);
                if (Math.Abs(fx) <= eps || Math.Abs(x - x2) <= eps) return x;

                if (fa * fx < 0f) { x2 = x; fb = fx; }
                else { x1 = x; fa = fx; }
            }
            return x;
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
            float x0 = 0.5f * (a + b);
            float fa = f(a);
            int n = 0;

            while (n++ < short.MaxValue)
            {
                float fx = f(x0);
                if (Math.Abs(fx) <= eps) return x0;

                float denom = (fa - fx);
                if (denom == 0f) break;

                float x1 = x0 - fx * (a - x0) / denom;
                if (Math.Abs(x1 - x0) <= eps) return x1;
                x0 = x1;
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
        private static Complex32 chord(IComplex32 f, Complex32 a, Complex32 b, float eps = 1e-8f)
        {
            Complex32 x0 = 0.5f * (a + b);
            Complex32 fa = f(a);
            int n = 0;

            while (n++ < short.MaxValue)
            {
                Complex32 fx = f(x0);
                if (Maths.Abs(fx) <= eps) return x0;

                Complex32 denom = (fa - fx);
                if (Maths.Abs(denom) == 0f) break;

                Complex32 x1 = x0 - fx * (a - x0) / denom;
                if (Maths.Abs(x1 - x0) <= eps) return x1;
                x0 = x1;
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
        private static Complex32 secan(IComplex32 f, Complex32 a, Complex32 b, float eps = 1e-8f)
        {
            Complex32 x1 = a, x2 = b;
            Complex32 f1 = f(x1), f2 = f(x2);
            int n = 0;

            while (Maths.Abs(f2) > eps && Maths.Abs(x2 - x1) > eps && n++ < short.MaxValue)
            {
                Complex32 denom = (f2 - f1);
                if (Maths.Abs(denom) == 0f) break;
                Complex32 x3 = x2 - (x2 - x1) * f2 / denom;

                x1 = x2; f1 = f2;
                x2 = x3; f2 = f(x2);
            }
            return x2;
        }
        #endregion
    }
}
