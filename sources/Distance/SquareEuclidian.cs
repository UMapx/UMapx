using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines square Euclidian distance.
    /// </summary>
    public class SquareEuclidian : DistanceBase, IDistance
    {
        #region Square-Euclidian distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0.0f;
            float u;

            for (int i = 0; i < n; i++)
            {
                u = p[i] - q[i];
                sum += u * u;
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
            float sum = 0.0f;
            float u;

            for (int i = 0; i < n; i++)
            {
                u = Maths.Abs(p[i] - q[i]);
                sum += u * u;
            }

            return sum;
        }
        #endregion
    }
}
