using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements an extremum search.
    /// <remarks>
    /// This class is a solution to the problem of finding the maximum and minimum points of the function F(x).
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Optimization
    {
        #region Private data
        private float eps;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements an extremum search.
        /// </summary>
        /// <param name="eps">Epsilon [0, 1]</param>
        public Optimization(float eps = 1e-8f)
        {
            this.Eps = eps;
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
        /// Returns the corresponding minimum of the function on the segment.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Start of line</param>
        /// <param name="b">End of line</param>
        /// <param name="max">Search maximum or minimum</param>
        /// <returns>Value</returns>
        public float Compute(IFloat function, float a, float b, bool max = false)
        {
            // max or min
            return (max) ? GoldenMax(function, a, b, this.eps) : GoldenMin(function, a, b, this.eps);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Golden-section search for the minimum of a unimodal function on [a, b].
        /// </summary>
        /// <remarks>
        /// - Assumes <paramref name="f"/> is continuous and unimodal on the closed interval [a, b].
        /// - Uses the golden ratio to reuse interval proportions for fast bracketing shrinkage.
        /// - Terminates when the bracket length |b - a| becomes smaller than <paramref name="eps"/> 
        ///   or when the iteration cap (short.MaxValue) is reached.
        /// - Returns the midpoint of the final bracket as the argmin approximation (not f at that point).
        /// 
        /// Note: This implementation recomputes f(x1) and f(x2) each iteration for clarity.
        /// It can be optimized to carry one evaluation forward per step.
        /// </remarks>
        /// <param name="f">Continuous objective function to minimize</param>
        /// <param name="a">Left endpoint of the search interval</param>
        /// <param name="b">Right endpoint of the search interval</param>
        /// <param name="eps">Absolute tolerance for the bracket length; stop when |b - a| &lt; eps.</param>
        /// <returns>
        /// Approximate minimizer x* ∈ [a, b] (the x-coordinate). To get the minimum value, evaluate f at the result.
        /// </returns>
        private static float GoldenMin(IFloat f, float a, float b, float eps = 1e-8f)
        {
            float x1, x2;

            for (int i = 0; i < short.MaxValue; i++)
            {
                x1 = b - (b - a) / Maths.Phi;
                x2 = a + (b - a) / Maths.Phi;

                if (f(x1) > f(x2))
                    a = x1;
                else
                    b = x2;
                if (Math.Abs(b - a) < eps)
                    break;
            }
            return (a + b) / 2;
        }
        /// <summary>
        /// Golden-section search for the maximum of a unimodal function on [a, b].
        /// </summary>
        /// <remarks>
        /// - Assumes <paramref name="f"/> is continuous and unimodal on [a, b].
        /// - Same scheme as <see cref="GoldenMin"/>, but with the inequality flipped to seek the maximum.
        /// - Terminates when |b - a| &lt; <paramref name="eps"/> or when the iteration cap is reached.
        /// - Returns the midpoint of the final bracket as the argmax approximation (not f at that point).
        /// </remarks>
        /// <param name="f">Continuous objective function to maximize</param>
        /// <param name="a">Left endpoint of the search interval</param>
        /// <param name="b">Right endpoint of the search interval</param>
        /// <param name="eps">Absolute tolerance for the bracket length; stop when |b - a| &lt; eps</param>
        /// <returns>
        /// Approximate maximizer x* ∈ [a, b] (the x-coordinate). To get the maximum value, evaluate f at the result.
        /// </returns>
        private static float GoldenMax(IFloat f, float a, float b, float eps = 1e-8f)
        {
            float x1, x2;

            for (int i = 0; i < short.MaxValue; i++)
            {
                x1 = b - (b - a) / Maths.Phi;
                x2 = a + (b - a) / Maths.Phi;

                if (f(x1) < f(x2))
                    a = x1;
                else
                    b = x2;
                if (Math.Abs(b - a) < eps)
                    break;
            }
            return (a + b) / 2;
        }
        #endregion
    }
}
