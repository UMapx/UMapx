using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Sokal-Michener distance.
    /// </summary>
    public class SokalMichener : DistanceBase, IDistance
    {
        #region Sokal-Michener distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int tf = 0;
            int ft = 0;
            int tt = 0;
            int ff = 0;
            int length = p.Length;

            for (int i = 0; i < length; i++)
            {
                if (p[i] == 1 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] == 1) ft++;
                if (p[i] == 1 && q[i] == 1) tt++;
                if (p[i] == 0 && q[i] == 0) ff++;
            }

            int n = tt + tf + ft + ff;
            return (tf + ft) / (float)n;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override Complex32 Compute(Complex32[] p, Complex32[] q)
        {
            int tf = 0;
            int ft = 0;
            int tt = 0;
            int ff = 0;
            int length = p.Length;

            for (int i = 0; i < length; i++)
            {
                if (p[i] == 1 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] == 1) ft++;
                if (p[i] == 1 && q[i] == 1) tt++;
                if (p[i] == 0 && q[i] == 0) ff++;
            }

            int n = tt + tf + ft + ff;
            return (tf + ft) / (float)n;
        }
        #endregion
    }
}
