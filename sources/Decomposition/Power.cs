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
        private double[] v;
        #endregion

        #region Power iteration components
        /// <summary>
        /// Initializes power iteration.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="iterations">Number of iterations</param>
        public Power(double[,] A, int iterations = 10)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("The matrix must be square");

            // eigenvalue power algorithm:
            int n = A.GetLength(0);
            this.v = Matrice.Rand(n);
            double[] w;
            double beta;

            // power iteration:
            for (int i = 0; i < iterations; i++)
            {
                // formula:
                // v[j] = (v[j-1] * A) / || v[j-1] * A ||
                w = Matrice.Dot(v, A);
                beta = Matrice.Norm(w);
                v = Matrice.Div(w, beta);
            }
            return;
        }
        /// <summary>
        /// Returns a vector of eigenvalues.
        /// </summary>
        public double[] V
        {
            get
            {
                return v;
            }
        }
        /// <summary>
        /// Returns the diagonalized matrix of eigenvalues.
        /// </summary>
        public double[,] J
        {
            get
            {
                return Matrice.Diag(v);
            }
        }
        #endregion
    }
}
