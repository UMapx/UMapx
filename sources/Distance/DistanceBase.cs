using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines the base distance class.
    /// </summary>
    public abstract class DistanceBase : IDistance
    {
        #region Distance methods
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public abstract float Compute(float[] p, float[] q);
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public abstract Complex32 Compute(Complex32[] p, Complex32[] q);
        /// <summary>
        /// Returns distance values.
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public virtual float[] Compute(float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = Compute(t, u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance values.
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public virtual Complex32[] Compute(Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = Compute(t, u);
            }

            return v;
        }
        #endregion
    }
}
