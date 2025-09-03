using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines RQ decomposition.
    /// <remarks>
    /// This is a matrix representation in the form of a product of two matrices: A = R * Q, 
    /// where Q is a unitary (or orthogonal) matrix, and R is an upper triangular matrix.
    /// RQ decomposition is one of the modifications of the QR algorithm.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/QR_decomposition
    /// </remarks>
    /// </summary>
    [Serializable]
    public class RQ
    {
        #region Private data
        private QR qr;
        private float[,] r;
        private float[,] q;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes RQ decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        public RQ(float[,] A)
        {
            qr = new QR(A.Flip(Direction.Vertical).Transponate());

            r = qr.R.Transponate();
            q = qr.Q.Transponate();

            r = r.Flip(Direction.Both);
            q = q.Flip(Direction.Vertical);
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Returns the lower triangular matrix R.
        /// </summary>
        public float[,] R
        {
            get
            {
                return this.r;
            }
        }
        /// <summary>
        /// Returns the orthogonal matrix Q.
        /// </summary>
        public float[,] Q
        {
            get
            {
                return this.q;
            }
        }
        #endregion
    }
}
