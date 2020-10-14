using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines LQ decomposition.
    /// <remarks>
    /// This is the representation of a matrix in the form of a product of two matrices: A = L * Q, where Q is a unitary (or orthogonal) matrix, and L is a lower triangular matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/QR_decomposition
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LQ
    {
        #region Private data
        private QR qr;
        private double[,] l;
        private double[,] q;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes LQ decomposition.
        /// </summary>
        /// <param name="A">Square matrix</param>
        public LQ(double[,] A)
        {
            qr = new QR(A.Transponate());

            l = qr.R.Transponate();
            q = qr.Q.Transponate();
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
