using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines decomposition with a cast to Hessenberg form.
    /// </summary>
    /// <remarks>
    /// This is a representation of a square matrix in the form of a product of three matrices: A = P * H * Pᵀ, 
    /// where H is the Hessenberg form and P is the unitary matrix.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hessenberg_matrix
    /// </remarks>
    [Serializable]
    public class Hessenberg
    {
        #region Private data
        private float[][] matrices;
        private float[][] hessenberg;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes decomposition to a Hessenberg form.
        /// </summary>
        /// <param name="A">Square matrix</param>
        public Hessenberg(float[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new ArgumentException("The matrix must be square");

            // Reduce to Hessenberg form.
            orthes(A);
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the unitary matrix.
        /// </summary>
        public float[,] P
        {
            get { return Jagged.FromJagged(matrices); }
        }
        /// <summary>
        /// Gets the Hessenberg form.
        /// </summary>
        public float[,] H
        {
            get { return Jagged.FromJagged(hessenberg); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Nonsymmetric reduction to Hessenberg form.
        /// This is derived from the Algol procedures orthes and ortran, by Martin and Wilkinson, 
        /// Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutines in EISPACK.
        /// </summary>
        /// <param name="A">Matrix</param>
        private void orthes(float[,] A)
        {
            // Properties
            int n = A.GetLength(0);
            this.matrices = Jagged.Zero(n, n);
            this.hessenberg = Jagged.ToJagged(A);
            float[] orthogonal = new float[n];

            int low = 0;
            int high = n - 1;
            int m, i, j;
            float scale, h, g, f;

            for (m = low + 1; m <= high - 1; m++)
            {
                // Scale column.

                scale = 0;
                for (i = m; i <= high; i++)
                    scale = scale + System.Math.Abs(hessenberg[i][m - 1]);

                if (scale != 0)
                {
                    // Compute Householder transformation.
                    h = 0;
                    for (i = high; i >= m; i--)
                    {
                        orthogonal[i] = hessenberg[i][m - 1] / scale;
                        h += orthogonal[i] * orthogonal[i];
                    }

                    g = (float)System.Math.Sqrt(h);
                    if (orthogonal[m] > 0) g = -g;

                    h = h - orthogonal[m] * g;
                    orthogonal[m] = orthogonal[m] - g;

                    // Apply Householder similarity transformation
                    // H = (I - u * u' / h) * H * (I - u * u') / h)
                    for (j = m; j < n; j++)
                    {
                        f = 0;
                        for (i = high; i >= m; i--)
                            f += orthogonal[i] * hessenberg[i][j];

                        f = f / h;
                        for (i = m; i <= high; i++)
                            hessenberg[i][j] -= f * orthogonal[i];
                    }

                    for (i = 0; i <= high; i++)
                    {
                        f = 0;
                        for (j = high; j >= m; j--)
                            f += orthogonal[j] * hessenberg[i][j];

                        f = f / h;
                        for (j = m; j <= high; j++)
                            hessenberg[i][j] -= f * orthogonal[j];
                    }

                    orthogonal[m] = scale * orthogonal[m];
                    hessenberg[m][m - 1] = scale * g;
                }
            }

            // Accumulate transformations (Algol's ortran).
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    matrices[i][j] = (i == j ? 1 : 0);

            for (m = high - 1; m >= low + 1; m--)
            {
                if (hessenberg[m][m - 1] != 0)
                {
                    for (i = m + 1; i <= high; i++)
                        orthogonal[i] = hessenberg[i][m - 1];

                    for (j = m; j <= high; j++)
                    {
                        g = 0;
                        for (i = m; i <= high; i++)
                            g += orthogonal[i] * matrices[i][j];

                        // float division avoids possible underflow.
                        g = (g / orthogonal[m]) / hessenberg[m][m - 1];
                        for (i = m; i <= high; i++)
                            matrices[i][j] += g * orthogonal[i];
                    }
                }
            }

            // final reduction:
            if (n > 2)
            {
                for (i = 0; i < n - 2; i++)
                {
                    for (j = i + 2; j < n; j++)
                    {
                        hessenberg[j][i] = 0;
                    }
                }
            }

        }
        #endregion
    }
}
