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
                throw new ArgumentException("The matrix must be square");

            // arnoldi decomposition
            ArnoldiDcmp(A, n, m);
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
        /// Computes Arnoldi transform for the input matrix.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="n">Dimension</param>
        /// <param name="m">Dimension</param>
        private void ArnoldiDcmp(float[,] a, int n, int m)
        {
            // vectors and matrices:
            this.q = new float[n, m];
            this.h = new float[n, m];
            float[,] p = new float[n, m + 1];
            float[] v, w;
            int i, j, k;

            // random 0-vector and norm:
            v = MatrixF.Rand(n);
            p = MatrixF.SetCol(p, v, 0);
            p = MatrixF.Div(p, MatrixF.Norm(v));

            // Start calculating
            // Arnoldi decomposition:
            for (k = 1; k <= m; k++)
            {
                // previous k-1-vector:
                v = MatrixF.Dot(MatrixF.GetCol(p, k - 1), a);

                for (j = 0; j < k; j++)
                {
                    // calculating α:
                    w = MatrixF.GetCol(p, j);
                    float alpha = MatrixF.Dot(w, v);
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
                    float beta = MatrixF.Norm(v);
                    p = MatrixF.SetCol(p, MatrixF.Div(v, beta), k);
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
