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
            int tt = 0;
            int tf = 0;
            int ft = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] != 0) tt++;
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
            }

            n = tt + tf + ft;
            return (n == 0) ? 0f : (tf + ft) / (float)n;
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
            int tt = 0;
            int tf = 0;
            int ft = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] != 0) tt++;
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
            }

            n = tt + tf + ft;
            return (n == 0) ? 0f : (tf + ft) / (float)n;
        }
        #endregion
    }
}
