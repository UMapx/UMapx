using System;
using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Chebyshev distance.
    /// </summary>
    public class Chebyshev : DistanceBase, IDistance
    {
        #region Chebyshev distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            float max = Math.Abs(p[0] - q[0]);
            float tmp;

            for (int k = 1; k < n; k++)
            {
                tmp = Math.Abs(p[k] - q[k]);
                max = tmp > max ? tmp : max;
            }

            return max;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override Complex32 Compute(Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            float max = Maths.Abs(p[0] - q[0]);
            float tmp;

            for (int k = 1; k < n; k++)
            {
                tmp = Maths.Abs(p[k] - q[k]);
                max = tmp > max ? tmp : max;
            }

            return max;
        }
        #endregion
    }
}
