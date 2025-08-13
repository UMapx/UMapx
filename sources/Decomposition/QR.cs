using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines QR decomposition.
    /// <remarks>
    /// This is a matrix representation in the form of a product of two matrices: A = Q * R, where Q is a unitary (or orthogonal) matrix, and R is an upper triangular matrix.
    /// QR decomposition is the basis of one of the search methods for eigenvectors and matrix numbers - the QR algorithm.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/QR_decomposition
    /// </remarks>
    /// </summary>
    [Serializable]
    public class QR
    {
        #region Private data
        private int m, n;
        private float[][] qr;
        private float[] diag;
        private float[,] q;
        private float[,] r;
        private float[,] h;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes QR decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        public QR(float[,] A)
        {
            qrdecomp(A);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Returns a matrix containing Householder reflection vectors.
        /// </summary>
        public float[,] H
        {
            get
            {
                return h;
            }
        }
        /// <summary>
        /// Returns the upper triangular matrix R.
        /// </summary>
        public float[,] R
        {
            get
            {
                return r;
            }
        }
        /// <summary>
        /// Returns the orthogonal matrix Q.
        /// </summary>
        public float[,] Q
        {
            get
            {
                return q;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        private void qrdecomp(float[,] A)
        {
            // params
            this.m = A.GetLength(0);
            this.n = A.GetLength(1);
            this.diag = new float[n];
            this.qr = Jagged.ToJagged(A);
            float nrm, s;
            int k, i, j;

            // Main loop.
            for (k = 0; k < n; k++)
            {
                // Compute 2-norm of k-th column without under/overflow.
                nrm = 0;

                for (i = k; i < m; i++)
                {
                    nrm = Maths.Hypotenuse(nrm, qr[i][k]);
                }

                if (nrm != 0.0)
                {
                    // Form k-th Householder vector.
                    if (qr[k][k] < 0)
                    {
                        nrm = -nrm;
                    }
                    for (i = k; i < m; i++)
                    {
                        qr[i][k] /= nrm; // Make v a unit vector
                    }
                    qr[k][k] += 1.0f; // + the (e)kth vector

                    // Apply transformation to remaining columns.
                    for (j = k + 1; j < n; j++) // For each column
                    {
                        s = 0.0f;
                        for (i = k; i < m; i++) // For each row
                        {
                            s += qr[i][k] * qr[i][j];
                        }

                        s = (-s) / qr[k][k]; // Unit vector product

                        for (i = k; i < m; i++) // For each row
                        {
                            qr[i][j] += s * qr[i][k];
                        }
                    }
                }

                diag[k] = -nrm;
            }

            // prepare Q matrix
            q = new float[m, n];

            for (k = n - 1; k >= 0; k--)
            {
                for (i = 0; i < m; i++)
                {
                    q[i, k] = 0.0f;
                }
                q[k, k] = 1.0f;
                for (j = k; j < n; j++)
                {
                    if (qr[k][k] != 0)
                    {
                        s = 0.0f;
                        for (i = k; i < m; i++)
                        {
                            s += qr[i][k] * q[i, j];
                        }
                        s = (-s) / qr[k][k];

                        for (i = k; i < m; i++)
                        {
                            q[i, j] += s * qr[i][k];
                        }
                    }
                }
            }

            // prepare R matrix
            r = new float[n, n];

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (i < j)
                    {
                        r[i, j] = qr[i][j];
                    }
                    else if (i == j)
                    {
                        r[i, j] = diag[i];
                    }
                    else
                    {
                        r[i, j] = 0.0f;
                    }
                }
            }

            // prepare H matrix
            h = new float[m, n];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (i >= j)
                    {
                        h[i, j] = qr[i][j];
                    }
                    else
                    {
                        h[i, j] = 0.0f;
                    }
                }
            }
        }
        #endregion
    }
}
