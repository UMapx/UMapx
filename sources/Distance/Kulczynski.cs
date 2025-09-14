using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Kulczynski distance.
    /// See https://en.wikipedia.org/wiki/Kulczynski_dissimilarity.
    /// </summary>
    public class Kulczynski : DistanceBase, IDistance
    {
        #region Kulczynski distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            float tf = 0;
            float ft = 0;
            float tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            return 0.5f * (tf / (tt + tf) + ft / (tt + ft));
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
            float tf = 0;
            float ft = 0;
            float tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            return 0.5f * (tf / (tt + tf) + ft / (tt + ft));
        }
        #endregion
    }
}
