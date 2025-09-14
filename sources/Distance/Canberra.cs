using System;
using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Canberra distance.
    /// </summary>
    public class Canberra : DistanceBase, IDistance
    {
        #region Canberra distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                float den = Math.Abs(p[i]) + Math.Abs(q[i]);
                if (den > 0)
                    sum += Math.Abs(p[i] - q[i]) / den;
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
            int n = p.Length;
            Complex32 sum = 0;

            for (int i = 0; i < n; i++)
            {
                float den = Maths.Abs(p[i]) + Maths.Abs(q[i]);
                if (den > 0)
                    sum += Maths.Abs(p[i] - q[i]) / den;
            }
            return sum;
        }
        #endregion
    }
}
