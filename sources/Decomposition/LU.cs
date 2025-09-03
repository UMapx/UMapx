using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines LU decomposition.
    /// <remarks>
    /// This is a representation of the square matrix A as the product of two matrices: A = L * U, where L is the lower triangular matrix, U is the upper triangular matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/LU_decomposition
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LU
    {
        #region Private data
        private float[][] lower;
        private float[][] upper;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes LU decomposition.
        /// </summary>
        /// <param name="A">Square matrix</param>
        public LU(float[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new ArgumentException("The matrix must be square");

            // LU-decomposition:
            LuDcmp(Jagged.ToJagged(A));
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the lower triangular matrix.
        /// </summary>
        public float[,] L
        {
            get { return Jagged.FromJagged(lower); }
        }
        /// <summary>
        /// Gets the upper triangular matrix.
        /// </summary>
        public float[,] U
        {
            get { return Jagged.FromJagged(upper); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Performs LU decomposition of the specified matrix.
        /// </summary>
        /// <param name="a">Square matrix to factorize</param>
        private void LuDcmp(float[][] a)
        {
            int i, j, k;
            int n = a.GetLength(0);
            float alpha, beta;
            this.upper = Jagged.Zero(n, n);
            this.lower = Jagged.Zero(n, n);

            for (i = 0; i < n; i++)
            {
                this.upper[i][i] = 1;
            }

            for (j = 0; j < n; j++)
            {
                for (i = j; i < n; i++)
                {
                    alpha = 0;
                    for (k = 0; k < j; k++)
                    {
                        alpha = alpha + this.lower[i][k] * this.upper[k][j];
                    }
                    this.lower[i][j] = a[i][j] - alpha;
                }

                beta = lower[j][j];

                for (i = j; i < n; i++)
                {
                    alpha = 0;
                    for (k = 0; k < j; k++)
                    {
                        alpha = alpha + this.lower[j][k] * this.upper[k][i];
                    }

                    if (beta != 0)
                    {
                        this.upper[j][i] = (a[j][i] - alpha) / beta;
                    }
                }
            }
        }
        #endregion
    }
}
