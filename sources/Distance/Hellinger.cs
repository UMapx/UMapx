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
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                float d = Maths.Sqrt(p[i]) - Maths.Sqrt(q[i]);
                sum += d * d;
            }

            return Maths.Sqrt(sum) / Maths.Sqrt2;
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
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                Complex32 diff = Maths.Sqrt(p[i]) - Maths.Sqrt(q[i]);
                float d = Maths.Abs(diff);
                sum += d * d;
            }

            return Maths.Sqrt(sum) / Maths.Sqrt2;
        }
        #endregion
    }
}
