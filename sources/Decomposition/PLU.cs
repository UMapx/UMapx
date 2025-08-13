using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines PLU decomposition.
    /// <remarks>
    /// This is the representation of a square matrix A as the product of three matrices: P * A = L * U,
    /// where P is the permutation matrix, L is the lower triangular matrix with unit diagonal, and U is
    /// the upper triangular matrix. The decomposition is computed with partial pivoting.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/LU_decomposition#LU_decomposition_with_partial_pivoting
    /// </remarks>
    /// </summary>
    [Serializable]
    public class PLU
    {
        #region Private data
        private float[][] perm;
        private float[][] lower;
        private float[][] upper;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes PLU decomposition.
        /// </summary>
        /// <param name="A">Square matrix</param>
        public PLU(float[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("The matrix must be square");

            // PLU-decomposition:
            pludecomp(Jagged.ToJagged(A));
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Gets the permutation matrix.
        /// </summary>
        public float[,] P
        {
            get { return Jagged.FromJagged(perm); }
        }
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
        /// Computes PLU decomposition with partial pivoting.
        /// </summary>
        /// <param name="a">Jagged array</param>
        private void pludecomp(float[][] a)
        {
            int n = a.GetLength(0);
            int i, j, k;
            lower = Jagged.Zero(n, n);
            upper = Jagged.Zero(n, n);
            perm = Jagged.Zero(n, n);

            for (i = 0; i < n; i++)
            {
                perm[i][i] = 1;
            }

            for (k = 0; k < n; k++)
            {
                // Find pivot
                int pivot = k;
                float max = Math.Abs(a[k][k]);
                for (i = k + 1; i < n; i++)
                {
                    float value = Math.Abs(a[i][k]);
                    if (value > max)
                    {
                        max = value;
                        pivot = i;
                    }
                }

                if (pivot != k)
                {
                    // Swap rows in a and permutation matrix
                    var tempRow = a[k]; a[k] = a[pivot]; a[pivot] = tempRow;
                    var tempPerm = perm[k]; perm[k] = perm[pivot]; perm[pivot] = tempPerm;

                    // Swap rows in lower matrix for previous columns
                    for (j = 0; j < k; j++)
                    {
                        float t = lower[k][j];
                        lower[k][j] = lower[pivot][j];
                        lower[pivot][j] = t;
                    }
                }

                // Compute U row
                for (j = k; j < n; j++)
                {
                    float sum = 0;
                    for (i = 0; i < k; i++)
                    {
                        sum += lower[k][i] * upper[i][j];
                    }
                    upper[k][j] = a[k][j] - sum;
                }

                // Compute L column
                lower[k][k] = 1;
                for (i = k + 1; i < n; i++)
                {
                    float sum = 0;
                    for (j = 0; j < k; j++)
                    {
                        sum += lower[i][j] * upper[j][k];
                    }
                    lower[i][k] = (a[i][k] - sum) / upper[k][k];
                }
            }
        }
        #endregion
    }
}

