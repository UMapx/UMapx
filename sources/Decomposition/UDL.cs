using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines UDL decomposition.
    /// <remarks>
    /// This is the representation of a symmetric square matrix as the product of three matrices: A = U * D * L, where U is the upper triangular matrix, D is the diagonal matrix, and L is the lower triangular matrix.
    /// This decomposition is a specific form of Cholesky decomposition.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class UDL
    {
        #region Private data
        private double[,] upper;
        private double[] diag;
        #endregion

        #region UDL components
        /// <summary>
        /// Initializes UDL decomposition.
        /// </summary>
        /// <param name="A">Square symmetric matrix</param>
        public UDL(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("The matrix must be square");

            udldecomp(A);
        }
        /// <summary>
        /// Returns the top triangular matrix.
        /// </summary>
        public double[,] U
        {
            get
            {
                return this.upper;
            }
        }
        /// <summary>
        /// Returns the diagonal matrix.
        /// </summary>
        public double[] D
        {
            get
            {
                return this.diag;
            }
        }
        /// <summary>
        /// Returns the lower triangular matrix.
        /// </summary>
        public double[,] L
        {
            get
            {
                return this.upper.Transponate();
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        private void udldecomp(double[,] a)
        {
            int i, j, k;
            int n = a.GetLength(0);
            this.upper = new double[n, n];
            this.diag = new double[n];
            double[][] p = Jagged.ToJagged(a);
            double alpha, beta, gamma;

            // Mathematics in science and engineering, v.128,
            // Factorization methods for discrete sequential estimation, Gerald J. Bierman.
            // UDU* factorization aglorithm.
            // 
            for (j = n - 1; j >= 1; j--)
            {
                gamma = p[j][j];
                diag[j] = gamma;
                alpha = 1.0 / gamma;

                for (k = 0; k < j; k++)
                {
                    beta = p[k][j];
                    upper[k, j] = alpha * beta;

                    for (i = 0; i <= k; i++)
                    {
                        p[i][k] -= beta * upper[i, j];
                    }
                }
            }
            diag[0] = p[0][0];

            // diagonal eyes:
            for (i = 0; i < n; i++)
            {
                upper[i, i] = 1.0;
            }
            return;
        }
        #endregion
    }
}
