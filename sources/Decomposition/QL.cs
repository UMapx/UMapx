using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines the QL decomposition of a square matrix.
    /// <remarks>
    /// This is a representation of a matrix in the form of a product of two matrices: A = Q * L, where Q is a unitary (or orthogonal) matrix and L is a lower triangular matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/QR_decomposition
    /// </remarks>
    /// </summary>
    [Serializable]
    public class QL
    {
        #region Private data
        private QR qr;
        private double[,] l;
        private double[,] q;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the QL decomposition of a square matrix.
        /// </summary>
        /// <param name="A">Square matrix</param>
        public QL(double[,] A)
        {
            qr = new QR(A.Flip(Direction.Horizontal));
            q = qr.Q.Flip(Direction.Horizontal);
            l = qr.R.Flip(Direction.Both);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Returns the lower triangular matrix L.
        /// </summary>
        public double[,] L
        {
            get
            {
                return this.l;
            }
        }
        /// <summary>
        /// Returns the orthogonal matrix Q.
        /// </summary>
        public double[,] Q
        {
            get
            {
                return this.q;
            }
        }
        #endregion
    }
}
