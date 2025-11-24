using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines non-negative matrix factorization.
    /// </summary>
    /// <remarks>
    /// This is a representation of a rectangular matrix A as the product of two matrices: A = W * H.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Non-negative_matrix_factorization
    /// </remarks>
    [Serializable]
    public class NMF
    {
        #region Private data
        private float[,] w;       // W is m x r (weights)
        private float[,] h;       // H is r x n (transformed data) (transposed)
        private readonly int n;   // number of input data vectors
        private readonly int m;   // dimension of input vector
        private readonly int r;   // dimension of output vector (reduced dimension)
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes non-negative matrix factorization.
        /// </summary>
        /// <param name="A">Non-negative matrix</param>
        /// <param name="r">The dimension of new matrices</param>
        /// <param name="iterations">Number of iterations</param>
        public NMF(float[,] A, int r, int iterations = 100)
        {
            this.m = A.GetLength(0);
            this.n = A.GetLength(1);
            this.r = r;

            // decompose
            NnmfDcmp(A, iterations);
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the left matrix.
        /// </summary>
        public float[,] W
        {
            get { return w; }
        }
        /// <summary>
        /// Gets the right matrix.
        /// </summary>
        public float[,] H
        {
            get { return h; }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Computes non-negative matrix factorization for input matrix.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="iterations">Iterations</param>
        private void NnmfDcmp(float[,] A, int iterations = 100)
        {
            // choose W and H randomly
            w = Matrice.Rand(m, r);
            h = Matrice.Rand(r, n);
            var Z = new float[r, r];

            // a small epsilon is added to the
            // denominator to avoid overflow.
            float eps = 10e-9f;
            int i, j, l, t;
            float s, d;

            // allocate
            var newW = new float[m, r];
            var newH = new float[r, n];

            for (t = 0; t < iterations; t++)
            {
                // Update H using the multiplicative
                // H = H .* (W'*A) ./ (W'*W*H + eps) 
                for (i = 0; i < r; i++)
                {
                    for (j = i; j < r; j++)
                    {
                        s = 0.0f;
                        for (l = 0; l < m; l++)
                            s += w[l, i] * w[l, j];
                        Z[i, j] = Z[j, i] = s;
                    }

                    for (j = 0; j < n; j++)
                    {
                        d = 0.0f;
                        for (l = 0; l < r; l++)
                            d += Z[i, l] * h[l, j];

                        s = 0.0f;
                        for (l = 0; l < m; l++)
                            s += w[l, i] * A[l, j];

                        newH[i, j] = h[i, j] * s / (d + eps);
                    }
                }

                // Update W using the multiplicative
                //   W = W .* (A*H') ./ (W*H*H' + eps)
                for (j = 0; j < r; j++)
                {
                    for (i = j; i < r; i++)
                    {
                        s = 0.0f;
                        for (l = 0; l < n; l++)
                            s += newH[i, l] * newH[j, l];
                        Z[i, j] = Z[j, i] = s;
                    }

                    for (i = 0; i < m; i++)
                    {
                        d = 0.0f;
                        for (l = 0; l < r; l++)
                            d += w[i, l] * Z[j, l];

                        s = 0.0f;
                        for (l = 0; l < n; l++)
                            s += A[i, l] * newH[j, l];

                        newW[i, j] = w[i, j] * s / (d + eps);
                    }
                }

                w = newW;
                h = newH;
            }
        }
        #endregion
    }
}
