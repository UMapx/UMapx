using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Euclidean distance.
    /// </summary>
    public class Euclidean : DistanceBase, IDistance
    {
        #region Euclidean distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            float sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Maths.Pow(p[k] - q[k], 2);
            }

            return Maths.Sqrt(sum);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override Complex32 Compute(Complex32[] p, Complex32[] q)
        {
            float sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                float d = Maths.Abs(p[k] - q[k]);
                sum += d * d;
            }

            return Maths.Sqrt(sum);
        }
        #endregion
    }
}
