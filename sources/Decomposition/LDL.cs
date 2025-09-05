using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines LDL decomposition.
    /// </summary>
    /// <remarks>
    /// This is a representation of a symmetric positive definite square matrix in the form of a product of three matrices: A = L * D * Lᵀ, 
    /// where L is a lower triangular matrix with strictly positive elements on the diagonal,
    /// and D is the diagonal matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition_2
    /// </remarks>
    [Serializable]
    public class LDL
    {
        #region Private data
        private Cholesky choldecomp;
        private Diagonal diagdecomp;
        private float[,] lower;
        private float[] diag;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes LDL decomposition.
        /// </summary>
        /// <param name="A">Square symmetric positive definite matrix</param>
        public LDL(float[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new ArgumentException("The matrix must be square");

            // LDL'-decomposition algorithm
            // Cholesky decomposition:
            choldecomp = new Cholesky(A);
            lower = choldecomp.L;

            // Diagonal decomposition:
            diagdecomp = new Diagonal(lower);
            lower = diagdecomp.B;
            diag = diagdecomp.D;

            // D = d^2:
            diag = Matrice.Mul(diag, diag);
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the lower triangular matrix L.
        /// </summary>
        public float[,] L
        {
            get { return lower; }
        }
        /// <summary>
        /// Gets the upper triangular matrix U.
        /// </summary>
        public float[,] U
        {
            get { return Matrice.Transpose(lower); }
        }
        /// <summary>
        /// Gets the diagonal matrix.
        /// </summary>
        public float[] D
        {
            get { return diag; }
        }
        #endregion
    }
}
