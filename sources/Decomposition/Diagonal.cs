using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines diagonal decomposition.
    /// <remarks>
    /// This is a representation of the square matrix A as the product of two matrices: A = B * D, where B is the Square matrix and D is the diagonal matrix.
    /// This decomposition is used to highlight diagonal matrices in other decompositions (for example, LDU-, LDL-decompositions).
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Diagonal
    {
        #region Private data
        private double[,] matrix;
        private double[] diag;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes diagonal decomposition.
        /// </summary>
        /// <param name="A">Square matrix</param>
        public Diagonal(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("The matrix must be square");

            int n = A.GetLength(0), i;
            this.diag = new double[n];

            for (i = 0; i < n; i++)
                diag[i] = A[i, i];

            this.matrix = Matrice.Dot(A, diag, true);
            return;
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Gets the square matrix.
        /// </summary>
        public double[,] B
        {
            get { return matrix; }
        }
        /// <summary>
        /// Gets the vector of diagonal elements.
        /// </summary>
        public double[] D
        {
            get { return diag; }
        }
        #endregion
    }
}
