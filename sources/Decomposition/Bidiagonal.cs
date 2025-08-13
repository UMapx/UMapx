using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines bidiagonal decomposition.
    /// <remarks>
    /// This is a representation of a matrix in the form A = U * B * V', where U and V are orthogonal matrices
    /// obtained using Householder transformations and B is a bidiagonal matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Bidiagonalization
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Bidiagonal
    {
        #region Private data
        private float[][] b;
        private float[][] u;
        private float[][] v;
        private int m, n;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes bidiagonal decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        public Bidiagonal(float[,] A)
        {
            decompose(A);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Gets the orthogonal matrix U.
        /// </summary>
        public float[,] U
        {
            get { return Jagged.FromJagged(u); }
        }
        /// <summary>
        /// Gets the bidiagonal matrix B.
        /// </summary>
        public float[,] B
        {
            get { return Jagged.FromJagged(b); }
        }
        /// <summary>
        /// Gets the orthogonal matrix V.
        /// </summary>
        public float[,] V
        {
            get { return Jagged.FromJagged(v); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Performs bidiagonal decomposition using Householder transformations.
        /// </summary>
        /// <param name="A">Matrix</param>
        private void decompose(float[,] A)
        {
            m = A.GetLength(0);
            n = A.GetLength(1);
            b = Jagged.ToJagged(A);
            u = Jagged.ToJagged(Matrice.Eye(m, m));
            v = Jagged.ToJagged(Matrice.Eye(n, n));
            int p = Math.Min(m, n);

            for (int k = 0; k < p; k++)
            {
                // -------- Left Householder (zero below diagonal in column k) --------
                float norm = 0f;
                for (int i = k; i < m; i++)
                    norm = Maths.Hypotenuse(norm, b[i][k]);

                if (norm != 0f)
                {
                    // Sign choice: make v0 = 1 + b[k][k]/norm well-conditioned
                    if (b[k][k] < 0f) norm = -norm;

                    // v := x/norm; v0 += 1  (store v in b[k..,k])
                    for (int i = k; i < m; i++)
                        b[i][k] /= norm;
                    b[k][k] += 1f;

                    float v0 = b[k][k];

                    // Apply to B on the left: B[:,j] += s * v, s = -(v^T B[:,j]) / v0
                    for (int j = k + 1; j < n; j++)
                    {
                        float s = 0f;
                        for (int i = k; i < m; i++)
                            s += b[i][k] * b[i][j];

                        s = -s / v0;

                        for (int i = k; i < m; i++)
                            b[i][j] += s * b[i][k];
                    }

                    // Accumulate U on the RIGHT: U = U * H_left  (update rows of U)
                    for (int i = 0; i < m; i++)
                    {
                        float s = 0f;
                        for (int t = k; t < m; t++)
                            s += u[i][t] * b[t][k];

                        s = -s / v0;

                        for (int t = k; t < m; t++)
                            u[i][t] += s * b[t][k];
                    }
                }

                // Save diagonal and zero below it in column k
                float dk = -norm;
                b[k][k] = dk;
                for (int i = k + 1; i < m; i++)
                    b[i][k] = 0f;

                if (k < n - 1)
                {
                    // -------- Right Householder (zero beyond superdiagonal in row k) --------
                    norm = 0f;
                    for (int j = k + 1; j < n; j++)
                        norm = Maths.Hypotenuse(norm, b[k][j]);

                    if (norm != 0f)
                    {
                        // Sign choice for row reflector
                        if (b[k][k + 1] < 0f) norm = -norm;

                        // w := row/norm; w0 += 1  (store w in b[k, k+1..])
                        for (int j = k + 1; j < n; j++)
                            b[k][j] /= norm;
                        b[k][k + 1] += 1f;

                        float w0 = b[k][k + 1];

                        // Apply to B on the right: B[i,:] += t * w^T, t = -(B[i,:]·w) / w0
                        for (int i = k + 1; i < m; i++)
                        {
                            float s = 0f;
                            for (int j = k + 1; j < n; j++)
                                s += b[i][j] * b[k][j];

                            s = -s / w0;

                            for (int j = k + 1; j < n; j++)
                                b[i][j] += s * b[k][j];
                        }

                        // Accumulate V on the RIGHT: V = V * H_right  (update rows of V)
                        for (int i = 0; i < n; i++)
                        {
                            float s = 0f;
                            for (int j = k + 1; j < n; j++)
                                s += v[i][j] * b[k][j];

                            s = -s / w0;

                            for (int j = k + 1; j < n; j++)
                                v[i][j] += s * b[k][j];
                        }
                    }

                    // Save superdiagonal and zero the rest to the right in row k
                    float ek = -norm;
                    b[k][k + 1] = ek;
                    for (int j = k + 2; j < n; j++)
                        b[k][j] = 0f;
                }
            }
        }
        #endregion
    }
}
