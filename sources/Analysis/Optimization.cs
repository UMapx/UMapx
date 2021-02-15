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
        /// <returns>float precision floating point number</returns>
        public float Compute(IFloat function, float a, float b, bool max = false)
        {
            // max or min
            return (max) ? goldenMax(function, a, b, this.eps) : goldenMin(function, a, b, this.eps);
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
        private static float goldenMin(IFloat f, float a, float b, float eps = 1e-8f)
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
        /// 
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        private static float goldenMax(IFloat f, float a, float b, float eps = 1e-8f)
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
