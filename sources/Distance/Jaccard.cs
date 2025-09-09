using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Jaccard distance.
    /// </summary>
    public class Jaccard : DistanceBase, IDistance
    {
        #region Jaccard distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            int inter = 0;
            int union = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 || q[i] != 0)
                {
                    if (Maths.Abs(p[i] - q[i]) < 1e-8f)
                        inter++;
                    union++;
                }
            }

            return (union == 0) ? 0 : 1.0f - (inter / (float)union);
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
            int inter = 0;
            int union = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 || q[i] != 0)
                {
                    if (Maths.Abs(p[i] - q[i]) < 1e-8f)
                        inter++;
                    union++;
                }
            }

            return (union == 0) ? 0 : 1.0f - (inter / (float)union);
        }
        #endregion
    }
}
