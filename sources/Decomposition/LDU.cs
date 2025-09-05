using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines LDU decomposition.
    /// </summary>
    /// <remarks>
    /// This is the representation of a square matrix A as the product of three matrices: A = L * D * U, 
    /// where L is the lower triangular matrix, D is the diagonal matrix, and U is the upper triangular matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/LU_decomposition
    /// </remarks>
    [Serializable]
    public class LDU
    {
        #region Private data
        private LU ludecomp;
        private Diagonal diagdecomp;
        private float[,] lower;
        private float[,] upper;
        private float[] diag;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes LDU decomposition.
        /// </summary>
        /// <param name="A">Square matrix</param>
        public LDU(float[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new ArgumentException("The matrix must be square");

            // LDU algorithm:
            // LU-decomposition:
            ludecomp = new LU(A);
            lower = ludecomp.L;
            upper = ludecomp.U;

            // Diagonal decomposition:
            diagdecomp = new Diagonal(lower);
            lower = diagdecomp.B;
            diag = diagdecomp.D;
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the lower triangular matrix.
        /// </summary>
        public float[,] L
        {
            get { return lower; }
        }
        /// <summary>
        /// Gets the upper triangular matrix.
        /// </summary>
        public float[,] U
        {
            get { return upper; }
        }
        /// <summary>
        /// Gets the vector of diagonal elements.
        /// </summary>
        public float[] D
        {
            get { return diag; }
        }
        #endregion
    }
}
