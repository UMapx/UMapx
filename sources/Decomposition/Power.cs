using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines power iteration.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Power_iteration
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Power
    {
        #region Private data
        private float[] v;
        #endregion

        #region Power iteration components
        /// <summary>
        /// Initializes power iteration.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="iterations">Number of iterations</param>
        public Power(float[,] A, int iterations = 10)
        {
            if (!MatrixF.IsSquare(A))
                throw new ArgumentException("The matrix must be square");

            // eigenvalue power algorithm:
            int n = A.GetLength(0);
            this.v = MatrixF.Rand(n);
            float[] w;
            float beta;

            // power iteration:
            for (int i = 0; i < iterations; i++)
            {
                // formula:
                // v[j] = (v[j-1] * A) / || v[j-1] * A ||
                w = MatrixF.Dot(v, A);
                beta = MatrixF.Norm(w);
                v = MatrixF.Div(w, beta);
            }
        }
        /// <summary>
        /// Returns a vector of eigenvalues.
        /// </summary>
        public float[] V
        {
            get
            {
                return v;
            }
        }
        /// <summary>
        /// Returns the diagonalized matrix of eigenvalues.
        /// </summary>
        public float[,] J
        {
            get
            {
                return MatrixF.Diag(v);
            }
        }
        #endregion
    }
}
