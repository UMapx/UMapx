using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Angular distance.
    /// </summary>
    public class Angular : DistanceBase, IDistance
    {
        #region Angular distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            float s = 0;
            float x = 0;
            float y = 0;

            for (int i = 0; i < n; i++)
            {
                s += p[i] * q[i];
                x += p[i] * p[i];
                y += q[i] * q[i];
            }

            float den = Maths.Sqrt(x) * Maths.Sqrt(y);
            if (den == 0)
                return Maths.Pi / 2;
            float cosine = s / den;
            // clamp cosine for numerical stability
            cosine = Maths.Max(-1f, Maths.Min(1f, cosine));

            return Maths.Acos(cosine);
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
            Complex32 s = 0;
            float x = 0;
            float y = 0;

            for (int i = 0; i < n; i++)
            {
                s += p[i] * q[i].Conjugate;
                float a = Maths.Abs(p[i]);
                float b = Maths.Abs(q[i]);
                x += a * a;
                y += b * b;
            }

            float den = Maths.Sqrt(x) * Maths.Sqrt(y);
            if (den == 0)
                return Maths.Pi / 2;
            float cosine = Maths.Abs(s) / den;
            // clamp cosine for numerical stability
            cosine = Maths.Max(-1f, Maths.Min(1f, cosine));

            return Maths.Acos(cosine);
        }
        #endregion
    }
}
