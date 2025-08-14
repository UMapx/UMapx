using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines QZ decomposition (generalized Schur decomposition).
    /// <remarks>
    /// This is a matrix representation for a pair of matrices A and B such that
    /// A = Q * S * Z' and B = Q * T * Z', where Q and Z are orthogonal and
    /// S and T are upper (quasi) triangular matrices.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/QZ_decomposition
    /// </remarks>
    /// </summary>
    [Serializable]
    public class QZ
    {
        #region Private data
        private readonly int n;
        private readonly float[,] q;
        private readonly float[,] z;
        private readonly float[,] s;
        private readonly float[,] t;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes QZ decomposition.
        /// </summary>
        /// <param name="A">Matrix A</param>
        /// <param name="B">Matrix B</param>
        /// <param name="eps">Epsilon [0, 1]</param>
        public QZ(float[,] A, float[,] B, float eps = 1e-16f)
        {
            if (A.GetLength(0) != A.GetLength(1))
                throw new ArgumentException("The matrix must be square");
            if (B.GetLength(0) != B.GetLength(1))
                throw new ArgumentException("The matrix must be square");
            if (A.GetLength(0) != B.GetLength(0) || A.GetLength(1) != B.GetLength(1))
                throw new ArgumentException("Matrices should be the same size");

            this.n = A.GetLength(0);

            var a = Jagged.ToJagged(A);
            var b = Jagged.ToJagged(B);
            var zMat = Jagged.Zero(n, n);
            for (int i = 0; i < n; i++) zMat[i][i] = 1f;
            int ierr = 0;

            GEVD.qzdecomp(a, b, Maths.Float(eps), zMat, ref ierr);

            if (ierr != 0)
                throw new Exception("QZ decomposition failed to converge");

            this.s = Jagged.FromJagged(a);
            this.t = Jagged.FromJagged(b);
            this.z = Jagged.FromJagged(zMat);

            float[,] bz = Matrice.Dot(B, this.z);
            float[,] tinv = InvertUpperTriangular(this.t);
            this.q = Matrice.Dot(bz, tinv);
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the orthogonal matrix Q.
        /// </summary>
        public float[,] Q { get { return q; } }
        /// <summary>
        /// Gets the orthogonal matrix Z.
        /// </summary>
        public float[,] Z { get { return z; } }
        /// <summary>
        /// Gets the quasi upper triangular matrix S.
        /// </summary>
        public float[,] S { get { return s; } }
        /// <summary>
        /// Gets the upper triangular matrix T.
        /// </summary>
        public float[,] T { get { return t; } }
        #endregion

        #region Private voids
        /// <summary>
        /// Inverts an upper triangular matrix.
        /// </summary>
        /// <param name="m">Matix</param>
        private static float[,] InvertUpperTriangular(float[,] m)
        {
            int n = m.GetLength(0);
            float[,] inv = new float[n, n];

            for (int i = n - 1; i >= 0; i--)
            {
                inv[i, i] = 1f / m[i, i];

                for (int j = i + 1; j < n; j++)
                {
                    float sum = 0f;
                    for (int k = i + 1; k <= j; k++)
                        sum += m[i, k] * inv[k, j];
                    inv[i, j] = -inv[i, i] * sum;
                }
            }

            return inv;
        }
        #endregion
    }
}
