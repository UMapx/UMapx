using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines Cholesky decomposition.
    /// </summary>
    /// <remarks>
    /// This is a representation of a symmetric positive definite square matrix in the form of a product: A = L * Lᵀ, 
    /// where L is a lower triangular matrix with strictly positive elements on the diagonal.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Cholesky_decomposition
    /// </remarks>
    [Serializable]
    public class Cholesky
    {
        #region Private data
        float[][] lower;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes Cholesky decomposition.
        /// </summary>
        /// <param name="A">Square symmetric positive definite matrix</param>
        public Cholesky(float[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new ArgumentException("The matrix must be square");

            // Cholesky decomposition:
            CholDcmp(Jagged.ToJagged(A));
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the lower triangular matrix L.
        /// </summary>
        public float[,] L
        {
            get { return Jagged.FromJagged(lower); }
        }
        /// <summary>
        /// Gets the upper triangular matrix U.
        /// </summary>
        public float[,] U
        {
            get { return Matrice.Transponate(L); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Computes Cholesky decomposition for the matrix A.
        /// </summary>
        /// <param name="a">Matrix</param>
        private void CholDcmp(float[][] a)
        {
            // Cholesky decomposition
            int n = a.GetLength(0);
            this.lower = new float[n][];
            float[] v, w, z, d = new float[n];
            float alpha;
            int j, i, k;

            // get diagonal elements
            for (i = 0; i < n; i++)
            {
                d[i] = a[i][i];
            }

            // do job
            for (j = 0; j < n; j++)
            {
                v = lower[j] = new float[n];
                z = a[j];

                for (i = 0; i <= j; i++)
                {
                    w = lower[i];
                    alpha = 0;

                    if (i == j)
                    {
                        for (k = 0; k < i; k++)
                        {
                            alpha += w[k] * w[k];
                        }

                        w[i] = (float)Math.Sqrt(d[i] - alpha);
                        lower[i] = w;
                    }
                    else
                    {
                        for (k = 0; k < i; k++)
                        {
                            alpha += w[k] * v[k];
                        }

                        v[i] = (z[i] - alpha) / w[i];
                    }
                }

                lower[j] = v;
            }
        }
        #endregion
    }
}
