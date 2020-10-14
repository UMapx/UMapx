using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines Lanczos transform.
    /// <remarks>
    /// This transformation is used to represent the symmetric matrix A as a product
    /// of three matrices: A = Q * T * Q', where T is a tridiagonal matrix, and Q is an orthogonal matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Lanczos_algorithm
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Lanczos
    {
        #region Private data
        private double[,] q;
        private double[,] t;
        #endregion

        #region Lanczos components
        /// <summary>
        /// Initializes Lanczos transformation.
        /// </summary>
        /// <param name="A">Symmetric matrix</param>
        /// <param name="full">Full reorthogonalization or not</param>
        public Lanczos(double[,] A, bool full = false)
        {
            // exception
            if (!Matrice.IsSymmetric(A))
                throw new Exception("The matrix must be symmetrical");

            // lanczos decomposition
            int n = A.GetLength(0);
            lanczos(A, n, full);
        }
        /// <summary>
        /// Returns the orthogonal matrix.
        /// </summary>
        public double[,] Q
        {
            get
            {
                return this.q;
            }
        }
        /// <summary>
        /// Returns a tridiagonal matrix.
        /// </summary>
        public double[,] T
        {
            get
            {
                return this.t;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="n"></param>
        /// <param name="full"></param>
        private void lanczos(double[,] a, int n, bool full)
        {
            // This function uses the Lanczos algorithm with full
            // re-orthogonalization to compute k x k symmetric tridiagonal
            // matrix T that approximates mat up to rank k with respect to
            // transformation Q. That is, A = Q * T * Q'.

            // params
            int i, j, y, k = n - 1;
            double[] v = Matrice.Rand(n), z;
            double[] u = v.Div(Matrice.Norm(v));
            this.q = Matrice.Eye(n, n).SetCol(u, 0);
            this.t = new double[n, n];
            double beta, alpha;

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

            return;
        }
        #endregion
    }
}
