using System;
using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Manhattan distance.
    /// </summary>
    public class Manhattan : DistanceBase, IDistance
    {
        #region Manhattan distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            float sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Math.Abs(p[k] - q[k]);
            }

            return sum;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override Complex32 Compute(Complex32[] p, Complex32[] q)
        {
            float sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Maths.Abs(p[k] - q[k]);
            }

            return sum;
        }
        #endregion
    }
}
