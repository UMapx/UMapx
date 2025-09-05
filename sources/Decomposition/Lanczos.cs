using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines Lanczos transform.
    /// </summary>
    /// <remarks>
    /// This transformation is used to represent the symmetric matrix A as a product
    /// of three matrices: A = Q * T * Qᵀ, where T is a tridiagonal matrix, and Q is an orthogonal matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Lanczos_algorithm
    /// </remarks>
    [Serializable]
    public class Lanczos
    {
        #region Private data
        private float[,] q;
        private float[,] t;
        #endregion

        #region Lanczos components
        /// <summary>
        /// Initializes Lanczos transformation.
        /// </summary>
        /// <param name="A">Symmetric matrix</param>
        /// <param name="full">Full reorthogonalization or not</param>
        public Lanczos(float[,] A, bool full = false)
        {
            // exception
            if (!Matrice.IsSymmetric(A))
                throw new ArgumentException("The matrix must be symmetrical");

            // lanczos decomposition
            int n = A.GetLength(0);
            LanczosDcmp(A, n, full);
        }
        /// <summary>
        /// Returns the orthogonal matrix.
        /// </summary>
        public float[,] Q
        {
            get
            {
                return this.q;
            }
        }
        /// <summary>
        /// Returns a tridiagonal matrix.
        /// </summary>
        public float[,] T
        {
            get
            {
                return this.t;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// This function uses the Lanczos algorithm with full
        /// re-orthogonalization to compute k x k symmetric tridiagonal
        /// matrix T that approximates mat up to rank k with respect to
        /// transformation Q. That is, A = Q * T * Q'.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="n">Dimension</param>
        /// <param name="full">Full or not</param>
        private void LanczosDcmp(float[,] a, int n, bool full)
        {
            // params
            int i, j, y, k = n - 1;
            float[] v = Matrice.Rand(n), z;
            float[] u = v.Div(Matrice.Norm(v));
            this.q = Matrice.Eye(n, n).SetCol(u, 0);
            this.t = new float[n, n];
            float beta, alpha;

            // do job
            for (i = 0; i < n; i++)
            {
                z = u.Dot(a);
                alpha = u.Dot(z);
                t[i, i] = alpha;

                // full re-orthogonalization (x1)
                for (j = 0; j <= i; j++)
                {
                    v = q.GetCol(j);
                    alpha = v.Dot(z);

                    for (y = 0; y < n; y++)
                        z[y] = z[y] - v[y] * alpha;
                }

                // full re-orthogonalization (x2)
                if (full)
                {
                    for (j = 0; j <= i; j++)
                    {
                        v = q.GetCol(j);
                        alpha = v.Dot(z);

                        for (y = 0; y < n; y++)
                            z[y] = z[y] - v[y] * alpha;
                    }
                }

                // tridiagonal matrix
                if (i < k)
                {
                    beta = Matrice.Norm(z);
                    u = z.Div(beta);
                    q = Matrice.SetCol(q, u, i + 1);
                    t[i, i + 1] = beta;
                    t[i + 1, i] = beta;
                }
            }
        }
        #endregion
    }
}
