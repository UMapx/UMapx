using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines QR decomposition with column pivoting (QRP).
    /// <remarks>
    /// This is the decomposition of a matrix A such that A * P = Q * R,
    /// where Q is an orthogonal matrix, R is an upper triangular matrix
    /// and P is a permutation matrix representing column pivots.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class QRP
    {
        #region Private data
        private int m, n;
        private float[][] qr;
        private float[] diag;
        private int[] piv;
        private float[,] q;
        private float[,] r;
        private float[,] h;
        private float[,] p;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes QRP decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        public QRP(float[,] A)
        {
            qrpdecomp(A);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Returns a matrix containing Householder reflection vectors.
        /// </summary>
        public float[,] H
        {
            get { return h; }
        }
        /// <summary>
        /// Returns the upper triangular matrix R.
        /// </summary>
        public float[,] R
        {
            get { return r; }
        }
        /// <summary>
        /// Returns the orthogonal matrix Q.
        /// </summary>
        public float[,] Q
        {
            get { return q; }
        }
        /// <summary>
        /// Returns the permutation matrix P.
        /// </summary>
        public float[,] P
        {
            get { return p; }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Computes QR decomposition with column pivoting.
        /// </summary>
        /// <param name="A">Matrix</param>
        private void qrpdecomp(float[,] A)
        {
            // params
            this.m = A.GetLength(0);
            this.n = A.GetLength(1);
            this.diag = new float[n];
            this.qr = Jagged.ToJagged(A);
            this.piv = new int[n];

            int i, j, k;
            float nrm, s;

            for (i = 0; i < n; i++)
            {
                piv[i] = i;
            }

            // Main loop with column pivoting
            for (k = 0; k < n; k++)
            {
                // Determine pivot column based on norms
                int pivCol = k;
                float maxNrm = 0.0f;
                for (j = k; j < n; j++)
                {
                    nrm = 0.0f;
                    for (i = k; i < m; i++)
                    {
                        nrm = Maths.Hypotenuse(nrm, qr[i][j]);
                    }
                    if (nrm > maxNrm)
                    {
                        maxNrm = nrm;
                        pivCol = j;
                    }
                }

                // Swap columns if needed
                if (pivCol != k)
                {
                    for (i = 0; i < m; i++)
                    {
                        float temp = qr[i][k];
                        qr[i][k] = qr[i][pivCol];
                        qr[i][pivCol] = temp;
                    }
                    int tp = piv[k];
                    piv[k] = piv[pivCol];
                    piv[pivCol] = tp;
                }

                // Compute 2-norm of k-th column
                nrm = 0.0f;
                for (i = k; i < m; i++)
                {
                    nrm = Maths.Hypotenuse(nrm, qr[i][k]);
                }

                if (nrm != 0.0f)
                {
                    // Form k-th Householder vector
                    if (qr[k][k] < 0)
                    {
                        nrm = -nrm;
                    }
                    for (i = k; i < m; i++)
                    {
                        qr[i][k] /= nrm;
                    }
                    qr[k][k] += 1.0f;

                    // Apply transformation to remaining columns
                    for (j = k + 1; j < n; j++)
                    {
                        s = 0.0f;
                        for (i = k; i < m; i++)
                        {
                            s += qr[i][k] * qr[i][j];
                        }
                        s = (-s) / qr[k][k];
                        for (i = k; i < m; i++)
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

            // prepare permutation matrix P
            p = new float[n, n];
            for (i = 0; i < n; i++)
            {
                p[piv[i], i] = 1.0f;
            }
        }
        #endregion
    }
}
