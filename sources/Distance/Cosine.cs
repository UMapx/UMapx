using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines cosine distance.
    /// </summary>
    public class Cosine : DistanceBase, IDistance
    {
        #region Contructor
        /// <summary>
        /// Initializes cosine distance.
        /// </summary>
        /// <param name="similarity">Use similarity formula or not</param>
        public Cosine(bool similarity = false) 
        {
            Similarity = similarity;
        }
        /// <summary>
        /// Use similarity formula or not.
        /// </summary>
        public bool Similarity { get; set; }
        #endregion

        #region Cosine distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Vector</param>
        /// <param name="b">Vector</param>
        /// <returns>Value; when either vector has zero magnitude, returns 0 if <see cref="Similarity"/> is true; otherwise returns 1</returns>
        public override float Compute(float[] p, float[] b)
        {
            int length = p.Length;
            float A = Matrice.Abs(p, false);
            float B = Matrice.Abs(b, false);
            float s = 0;

            for (int i = 0; i < length; i++)
                s += p[i] * b[i];

            if (A == 0 || B == 0)
                return Similarity ? 0f : 1f;

            var similarity = s / (A * B);
            return Similarity ? similarity : 1f - similarity;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Vector</param>
        /// <param name="b">Vector</param>
        /// <returns>Value; when either vector has zero magnitude, returns 0 if <see cref="Similarity"/> is true; otherwise returns 1</returns>
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
                return Similarity ? 0f : 1f;

            float similarity = Maths.Abs(s) / den;
            return Similarity ? similarity : 1f - similarity;
        }
        #endregion
    }
}
