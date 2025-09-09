using System;
using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Bray-Curtis distance.
    /// </summary>
    public class BrayCurtis : DistanceBase, IDistance
    {
        #region Bray-Curtis distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            float x = 0;
            float y = 0;

            for (int i = 0; i < n; i++)
            {
                y += Math.Abs(p[i] - q[i]);
                x += Math.Abs(p[i] + q[i]);
            }

            return y / x;
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
            Complex32 x = 0;
            Complex32 y = 0;

            for (int i = 0; i < n; i++)
            {
                y += Maths.Abs(p[i] - q[i]);
                x += Maths.Abs(p[i] + q[i]);
            }

            return y / x;
        }
        #endregion
    }
}
