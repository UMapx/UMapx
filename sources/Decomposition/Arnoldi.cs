using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines Arnoldi transform.
    /// <remarks>
    /// This transformation is used to reduce the square matrix to the Hessenberg form.
    /// The matrix A is represented as the product of three matrices: A = Q * H * Q', where H is the upper Hessenberg triangular matrix, Q is the orthogonal matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Arnoldi_iteration
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Arnoldi
    {
        #region Private data
        private float[,] q;
        private float[,] h;
        #endregion

        #region Arnoldi components
        /// <summary>
        /// Initializes Arnoldi transformation.
        /// </summary>
        /// <param name="A">Square matrix</param>
        public Arnoldi(float[,] A)
        {
            // matrix properties
            int n = A.GetLength(0);
            int m = A.GetLength(1);

            if (n != m)
                throw new Exception("The matrix must be square");

            // arnoldi decomposition
            arnoldi(A, n, m);
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
        /// Returns the upper triangular Hessenberg matrix.
        /// </summary>
        public float[,] H
        {
            get
            {
                return this.h;
            }
        }
        #endregion

        #region Private data
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="n"></param>
        /// <param name="m"></param>
        private void arnoldi(float[,] a, int n, int m)
        {
            // vectors and matrices:
            this.q = new float[n, m];
            this.h = new float[n, m];
            float[,] p = new float[n, m + 1];
            float[] v, w;
            float alpha = 0, beta = 0;
            int i, j, k;

            // random 0-vector and norm:
            v = Matrice.Rand(n);
            p = Matrice.SetCol(p, v, 0);
            p = Matrice.Div(p, Matrice.Norm(v));

            // Start calculating
            // Arnoldi decomposition:
            for (k = 1; k <= m; k++)
            {
                // previous k-1-vector:
                v = Matrice.Dot(Matrice.GetCol(p, k - 1), a);

                for (j = 0; j < k; j++)
                {
                    // calculating α:
                    w = Matrice.GetCol(p, j);
                    alpha = Matrice.Dot(w, v);
                    h[j, k - 1] = alpha;

                    // finding k-vector:
                    for (i = 0; i < n; i++)
                    {
                        v[i] -= w[i] * alpha;
                    }
                }

                // transform:
                if (k < m)
                {
                    // calculating β:
                    beta = Matrice.Norm(v);
                    p = Matrice.SetCol(p, Matrice.Div(v, beta), k);
                    h[k, k - 1] = beta;
                }
            }

            // result Q-matrix:
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < m; j++)
                {
                    q[i, j] = p[i, j];
                }
            }
        }
        #endregion
    }
}
