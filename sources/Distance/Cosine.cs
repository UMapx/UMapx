using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines cosine distance.
    /// </summary>
    public class Cosine : DistanceBase, IDistance
    {
        #region Cosine distance
        /// <summary>
        /// Returns similarity function of two vectors.
        /// </summary>
        /// <param name="p">Vector</param>
        /// <param name="b">Vector</param>
        /// <returns>Value; returns 0 if either vector has zero magnitude</returns>
        public override float Compute(float[] p, float[] b)
        {
            int length = p.Length;
            float A = Matrice.Abs(p, false);
            float B = Matrice.Abs(b, false);
            float s = 0;

            for (int i = 0; i < length; i++)
                s += p[i] * b[i];

            if (A == 0 || B == 0)
                return 0;

            return s / (A * B);
        }
        /// <summary>
        /// Returns similarity function of two vectors.
        /// </summary>
        /// <param name="p">Vector</param>
        /// <param name="b">Vector</param>
        /// <returns>Value; returns 0 if either vector has zero magnitude</returns>
        public override Complex32 Compute(Complex32[] p, Complex32[] b)
        {
            int length = p.Length;
            Complex32 s = 0;
            float A = 0;
            float B = 0;

            for (int i = 0; i < length; i++)
            {
                s += p[i] * b[i].Conjugate;
                float ap = Maths.Abs(p[i]);
                float bp = Maths.Abs(b[i]);
                A += ap * ap;
                B += bp * bp;
            }

            float den = Maths.Sqrt(A) * Maths.Sqrt(B);
            if (den == 0)
                return 0;

            return s / den;
        }
        #endregion
    }
}
