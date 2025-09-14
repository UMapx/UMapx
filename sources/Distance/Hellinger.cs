using System;
using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Hellinger distance.
    /// </summary>
    public class Hellinger : DistanceBase, IDistance
    {
        #region Hellinger distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array of non-negative values</param>
        /// <param name="q">Array of non-negative values</param>
        /// <returns>Value</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] < 0 || q[i] < 0)
                    throw new ArgumentException("Arrays must contain non-negative values");

                float d = Maths.Sqrt(p[i]) - Maths.Sqrt(q[i]);
                sum += d * d;
            }

            return Maths.Sqrt(sum) / Maths.Sqrt2;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array of complex numbers with non-negative real parts</param>
        /// <param name="q">Array of complex numbers with non-negative real parts</param>
        /// <returns>Value</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        public override Complex32 Compute(Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i].Real < 0 || q[i].Real < 0)
                    throw new ArgumentException("Arrays must contain non-negative real parts");

                Complex32 diff = Maths.Sqrt(p[i]) - Maths.Sqrt(q[i]);
                float d = Maths.Abs(diff);
                sum += d * d;
            }

            return Maths.Sqrt(sum) / Maths.Sqrt2;
        }
        #endregion
    }
}
