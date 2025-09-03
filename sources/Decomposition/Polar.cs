using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines polar decomposition.
    /// <remarks>
    /// This is a representation of a rectangular matrix A in the form of a product of two matrices: A = U * P, 
    /// where U is a unitary matrix, P is a positive definite matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Polar_decomposition
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Polar
    {
        #region Private data
        private readonly SVD svd;
        readonly float[,] u;
        readonly float[,] p;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes polar decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="iterations">Number of iterations</param>
        public Polar(float[,] A, int iterations = 10)
        {
            svd = new SVD(A, iterations);

            float[,] U = svd.U;
            float[,] V = svd.V;
            float[,] H = V.Transponate();
            float[] S = svd.S;

            u = U.Dot(H); p = V.Dot(S).Dot(H);
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the unitary matrix.
        /// </summary>
        public float[,] U
        {
            get
            {
                return this.u;
            }
        }
        /// <summary>
        /// Gets a positive definite matrix.
        /// </summary>
        public float[,] P
        {
            get
            {
                return this.p;
            }
        }
        #endregion
    }
}
