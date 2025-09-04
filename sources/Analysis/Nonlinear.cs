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
                this.eps = MathsF.Float(value);
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
                    return Nonlinear.Chord(function, a, b, this.eps);
                case NonlinearMethod.FalsePosition:
                    return Nonlinear.Falpo(function, a, b, this.eps);
                case NonlinearMethod.Bisection:
                    return Nonlinear.Bisec(function, a, b, this.eps);
                default:
                    return Nonlinear.Secan(function, a, b, this.eps);
            }
        }
        /// <summary>
        /// Gets the root value of a nonlinear equation.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Start of line</param>
        /// <param name="b">End of line</param>
        /// <returns>Value</returns>
        public ComplexF Compute(IComplexF function, ComplexF a, ComplexF b)
        {
            // chose method of nonlinear
            switch (method)
            {
                case NonlinearMethod.Chord:
                    return Nonlinear.Chord(function, a, b, this.eps);
                case NonlinearMethod.FalsePosition:
                    throw new NotSupportedException("False position is not defined for complex-valued functions");
                case NonlinearMethod.Bisection:
                    throw new NotSupportedException("Bisection is not defined for complex-valued functions");
                default:
                    return Nonlinear.Secan(function, a, b, this.eps);
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Bisection method for solving f(x) = 0 on a bracketing interval [a, b].
        /// </summary>
        /// <remarks>
        /// - Requires a sign change on the endpoints: f(a) * f(b) ≤ 0 (throws if violated).<br/>
        /// - Guarantees convergence to a root in [a, b] for continuous f (linear convergence).<br/>
        /// - Stops when the interval width or |f(mid)| is below <paramref name="eps"/>, or on iteration cap.
        /// - Returns the final midpoint or a short secant refinement if available.
        /// </remarks>
        /// <param name="f">Scalar continuous function</param>
        /// <param name="a">Left endpoint of the initial bracket</param>
        /// <param name="b">Right endpoint of the initial bracket</param>
        /// <param name="eps">Absolute tolerance for both x-interval and residual checks</param>
        /// <returns>Approximate root in [a, b]</returns>
        /// <exception cref="ArgumentException">If f(a) and f(b) have the same strict sign</exception>
        private static float Bisec(IFloat f, float a, float b, float eps = 1e-8f)
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
        /// Secant method (open method) for solving f(x) = 0 using two initial guesses.
        /// </summary>
        /// <remarks>
        /// - Does not require a bracket; may fail if the sequence diverges or hits a zero derivative surrogate.<br/>
        /// - Superlinear convergence near a simple root; sensitive to starting points.<br/>
        /// - Stops on small residual |f(x)|, small step |Δx|, or iteration cap; returns the last iterate.
        /// </remarks>
        /// <param name="f">Scalar continuous function</param>
        /// <param name="a">First initial guess</param>
        /// <param name="b">Second initial guess</param>
        /// <param name="eps">Absolute tolerance for residual and step size</param>
        /// <returns>Approximate root</returns>
        private static float Secan(IFloat f, float a, float b, float eps = 1e-8f)
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
        /// Regula Falsi (false position) method for solving f(x) = 0 on a bracketing interval [a, b].
        /// </summary>
        /// <remarks>
        /// - Maintains a bracket with a sign change at every step (robust).<br/>
        /// - Can stagnate if one endpoint changes very slowly; convergence is at least linear.<br/>
        /// - Stops on small residual |f(x)|, small step |Δx|, or iteration cap; returns the last secant point.
        /// </remarks>
        /// <param name="f">Scalar continuous function</param>
        /// <param name="a">Left endpoint of the initial bracket</param>
        /// <param name="b">Right endpoint of the initial bracket</param>
        /// <param name="eps">Absolute tolerance for residual and step size</param>
        /// <returns>Approximate root in [a, b]</returns>
        /// <exception cref="ArgumentException">If f(a) and f(b) have the same strict sign</exception>
        private static float Falpo(IFloat f, float a, float b, float eps = 1e-8f)
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
        /// One-endpoint chord method (fixed endpoint secant) for solving f(x) = 0.
        /// </summary>
        /// <remarks>
        /// - Uses fixed endpoint at <paramref name="a"/> and forms a secant with the current iterate.<br/>
        /// - Does not enforce bracketing; convergence is problem-dependent and may be slow or fail.<br/>
        /// - Stops on small residual |f(x)|, small step |Δx|, or iteration cap.
        /// </remarks>
        /// <param name="f">Scalar continuous function</param>
        /// <param name="a">Fixed endpoint used in each chord</param>
        /// <param name="b">Second endpoint used only to initialize the first iterate</param>
        /// <param name="eps">Absolute tolerance for residual and step size</param>
        /// <returns>Approximate root</returns>
        private static float Chord(IFloat f, float a, float b, float eps = 1e-8f)
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
        /// One-endpoint chord method (fixed endpoint secant) for solving f(x) = 0.
        /// </summary>
        /// <remarks>
        /// - Uses fixed endpoint at <paramref name="a"/> and forms a secant with the current iterate.<br/>
        /// - Does not enforce bracketing; convergence is problem-dependent and may be slow or fail.<br/>
        /// - Stops on small residual |f(x)|, small step |Δx|, or iteration cap.
        /// </remarks>
        /// <param name="f">Scalar continuous function</param>
        /// <param name="a">Fixed endpoint used in each chord</param>
        /// <param name="b">Second endpoint used only to initialize the first iterate</param>
        /// <param name="eps">Absolute tolerance for residual and step size</param>
        /// <returns>Approximate root</returns>
        private static ComplexF Chord(IComplexF f, ComplexF a, ComplexF b, float eps = 1e-8f)
        {
            ComplexF x0 = 0.5f * (a + b);
            ComplexF fa = f(a);
            int n = 0;

            while (n++ < short.MaxValue)
            {
                ComplexF fx = f(x0);
                if (MathsF.Abs(fx) <= eps) return x0;

                ComplexF denom = (fa - fx);
                if (MathsF.Abs(denom) == 0f) break;

                ComplexF x1 = x0 - fx * (a - x0) / denom;
                if (MathsF.Abs(x1 - x0) <= eps) return x1;
                x0 = x1;
            }
            return x0;
        }
        /// <summary>
        /// Secant method (open method) for solving f(x) = 0 using two initial guesses.
        /// </summary>
        /// <remarks>
        /// - Does not require a bracket; may fail if the sequence diverges or hits a zero derivative surrogate.<br/>
        /// - Superlinear convergence near a simple root; sensitive to starting points.<br/>
        /// - Stops on small residual |f(x)|, small step |Δx|, or iteration cap; returns the last iterate.
        /// </remarks>
        /// <param name="f">Scalar continuous function</param>
        /// <param name="a">First initial guess</param>
        /// <param name="b">Second initial guess</param>
        /// <param name="eps">Absolute tolerance for residual and step size</param>
        /// <returns>Approximate root</returns>
        private static ComplexF Secan(IComplexF f, ComplexF a, ComplexF b, float eps = 1e-8f)
        {
            ComplexF x1 = a, x2 = b;
            ComplexF f1 = f(x1), f2 = f(x2);
            int n = 0;

            while (MathsF.Abs(f2) > eps && MathsF.Abs(x2 - x1) > eps && n++ < short.MaxValue)
            {
                ComplexF denom = (f2 - f1);
                if (MathsF.Abs(denom) == 0f) break;
                ComplexF x3 = x2 - (x2 - x1) * f2 / denom;

                x1 = x2; f1 = f2;
                x2 = x3; f2 = f(x2);
            }
            return x2;
        }
        #endregion
    }
}
