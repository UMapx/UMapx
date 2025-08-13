using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines bidiagonal decomposition.
    /// <remarks>
    /// This is a representation of a matrix in the form A = U * B * V', where U and V are orthogonal matrices
    /// obtained using Householder transformations and B is a bidiagonal matrix.
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
                // Left Householder transformation
                float norm = 0;
                for (int i = k; i < m; i++)
                    norm = Maths.Hypotenuse(norm, b[i][k]);

                if (norm != 0)
                {
                    if (b[k][k] > 0) norm = -norm;
                    for (int i = k; i < m; i++) b[i][k] /= norm;
                    b[k][k] += 1f;

                    for (int j = k + 1; j < n; j++)
                    {
                        float s = 0;
                        for (int i = k; i < m; i++) s += b[i][k] * b[i][j];
                        s = -s / b[k][k];
                        for (int i = k; i < m; i++) b[i][j] += s * b[i][k];
                    }

                    for (int j = 0; j < m; j++)
                    {
                        float s = 0;
                        for (int i = k; i < m; i++) s += b[i][k] * u[i][j];
                        s = -s / b[k][k];
                        for (int i = k; i < m; i++) u[i][j] += s * b[i][k];
                    }
                }
                float dk = -norm;
                b[k][k] = dk;
                for (int i = k + 1; i < m; i++) b[i][k] = 0f;

                if (k < n - 1)
                {
                    // Right Householder transformation
                    norm = 0;
                    for (int j = k + 1; j < n; j++)
                        norm = Maths.Hypotenuse(norm, b[k][j]);

                    if (norm != 0)
                    {
                        if (b[k][k + 1] > 0) norm = -norm;
                        for (int j = k + 1; j < n; j++) b[k][j] /= norm;
                        b[k][k + 1] += 1f;

                        for (int i = k + 1; i < m; i++)
                        {
                            float s = 0;
                            for (int j = k + 1; j < n; j++) s += b[i][j] * b[k][j];
                            s = -s / b[k][k + 1];
                            for (int j = k + 1; j < n; j++) b[i][j] += s * b[k][j];
                        }

                        for (int j = 0; j < n; j++)
                        {
                            float s = 0;
                            for (int i = k + 1; i < n; i++) s += b[k][i] * v[i][j];
                            s = -s / b[k][k + 1];
                            for (int i = k + 1; i < n; i++) v[i][j] += s * b[k][i];
                        }
                    }
                    float ek = -norm;
                    b[k][k + 1] = ek;
                    for (int j = k + 2; j < n; j++) b[k][j] = 0f;
                }
            }
        }
        #endregion
    }
}
