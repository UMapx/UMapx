using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines Householder transformation.
    /// </summary>
    /// <remarks>
    /// This is a linear transformation H (u) of the vector space V, which describes its mapping with respect to the hyperplane,
    /// which passes through the origin. It was proposed in 1958 by the American mathematician Elston Scott Householder. 
    /// Widely used in linear algebra for QR decomposition of a matrix.
    /// In addition, the Householder transform is actively used for orthogonalization of bases; ultimately, the Householder matrix has the following properties:
    /// H = Hᵀ, Hᵀ * H = I; det(H) = -1.
    /// In this class, two types of the Householder transform are implemented: reduction to a three-diagonal matrix and construction of the 
    /// Householder matrix from a given vector.
    /// In the first case, the original square matrix is defined as: A = H * T * Hᵀ.
    /// More information can be found on the website: 
    /// https://en.wikipedia.org/wiki/Householder_transformation
    /// </remarks>
    [Serializable]
    public class Householder
    {
        #region Private data
        private int n;
        private float[] Re, Im;
        private float[][] matrices;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes Householder transformation.
        /// </summary>
        /// <param name="v">Array</param>
        public Householder(float[] v)
        {
            // properties:
            this.n = v.Length;
            this.Re = Matrice.One(n);
            this.Im = new float[n];

            // reflection to 
            // Householder matrix:
            hmatx(v);
        }
        /// <summary>
        /// Initializes Householder transformation.
        /// </summary>
        /// <param name="A">Symmetric matrix</param>
        public Householder(float[,] A)
        {
            if (!Matrice.IsSymmetric(A))
                throw new ArgumentException("The matrix must be symmetric");

            // properties:
            this.n = A.GetLength(0);
            this.Re = new float[n];
            this.Im = new float[n];

            // reduction to 
            // tridiagonalization matrix:
            tred2(A);
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Returns the Householder matrix.
        /// </summary>
        public float[,] H
        {
            get
            {
                return Jagged.FromJagged(matrices);
            }
        }
        /// <summary>
        /// Gets the diagonal matrix.
        /// </summary>
        public float[,] T
        {
            get
            {
                float[,] D = new float[n, n];
                int i;

                // diagonal:
                for (i = 0; i < n; i++)
                {
                    D[i, i] = Re[i];
                }
                // diagonal left and right 
                // sides:
                for (i = 1; i < n; i++)
                {
                    D[i - 1, i] = Im[i];
                    D[i, i - 1] = Im[i];
                }

                return D;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Builds a full Householder reflector matrix for the input vector.
        /// </summary>
        /// <param name="v">Vector</param>
        private void hmatx(float[] v)
        {
            // [1] Alston S. Householder, "Unitary Triangularization of a Nonsymmetric Matrix", 
            // Journal of the ACM 5, 339-242, 1958;
            // [2] G. W. Stewart, Matrix Algorithms: Volume 1: Basic Decompositions, SIAM, xix+458, 1998.
            // 
            // Get Householder vector:
            float[] w = Matrice.Householder(v);

            // Get Householder matrix:
            int n = w.Length, i, j;
            float[] z;
            this.matrices = new float[n][];

            // M = I - w * w':
            for (i = 0; i < n; i++)
            {
                // eye vector
                z = new float[n];
                z[i] = 1.0f;

                for (j = 0; j < n; j++)
                {
                    z[j] -= w[i] * w[j];
                }

                matrices[i] = z;
            }
        }
        /// <summary>
        /// Symmetric Householder reduction to tridiagonal form.
        /// This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch, and Wilkinson, 
        /// Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
        /// </summary>
        /// <param name="a">Matrix</param>
        private void tred2(float[,] a)
        {
            int i, j, k;
            this.matrices = Jagged.ToJagged(a);

            for (j = 0; j < n; j++)
            {
                Re[j] = matrices[n - 1][j];
            }

            float scale, h, f, g, hh;

            // Householder reduction to tridiagonal form.
            for (i = n - 1; i > 0; i--)
            {
                // Scale to avoid under/overflow.
                scale = 0;
                h = 0;
                for (k = 0; k < i; k++)
                    scale = scale + System.Math.Abs(Re[k]);

                if (scale == 0)
                {
                    Im[i] = Re[i - 1];
                    for (j = 0; j < i; j++)
                    {
                        Re[j] = matrices[i - 1][j];
                        matrices[i][j] = 0;
                        matrices[j][i] = 0;
                    }
                }
                else
                {
                    // Generate Householder Matrice.
                    for (k = 0; k < i; k++)
                    {
                        Re[k] /= scale;
                        h += Re[k] * Re[k];
                    }

                    f = Re[i - 1];
                    g = (float)System.Math.Sqrt(h);
                    if (f > 0) g = -g;

                    Im[i] = scale * g;
                    h = h - f * g;
                    Re[i - 1] = f - g;
                    for (j = 0; j < i; j++)
                        Im[j] = 0;

                    // Apply similarity transformation to remaining columns.
                    for (j = 0; j < i; j++)
                    {
                        f = Re[j];
                        matrices[j][i] = f;
                        g = Im[j] + matrices[j][j] * f;
                        for (k = j + 1; k <= i - 1; k++)
                        {
                            g += matrices[k][j] * Re[k];
                            Im[k] += matrices[k][j] * f;
                        }
                        Im[j] = g;
                    }

                    f = 0;
                    for (j = 0; j < i; j++)
                    {
                        Im[j] /= h;
                        f += Im[j] * Re[j];
                    }

                    hh = f / (h + h);
                    for (j = 0; j < i; j++)
                        Im[j] -= hh * Re[j];

                    for (j = 0; j < i; j++)
                    {
                        f = Re[j];
                        g = Im[j];
                        for (k = j; k <= i - 1; k++)
                            matrices[k][j] -= (f * Im[k] + g * Re[k]);

                        Re[j] = matrices[i - 1][j];
                        matrices[i][j] = 0;
                    }
                }
                Re[i] = h;
            }

            // Accumulate transformations.
            for (i = 0; i < n - 1; i++)
            {
                matrices[n - 1][i] = matrices[i][i];
                matrices[i][i] = 1;
                h = Re[i + 1];
                if (h != 0)
                {
                    for (k = 0; k <= i; k++)
                        Re[k] = matrices[k][i + 1] / h;

                    for (j = 0; j <= i; j++)
                    {
                        g = 0;
                        for (k = 0; k <= i; k++)
                            g += matrices[k][i + 1] * matrices[k][j];
                        for (k = 0; k <= i; k++)
                            matrices[k][j] -= g * Re[k];
                    }
                }

                for (k = 0; k <= i; k++)
                    matrices[k][i + 1] = 0;
            }

            for (j = 0; j < n; j++)
            {
                Re[j] = matrices[n - 1][j];
                matrices[n - 1][j] = 0;
            }

            matrices[n - 1][n - 1] = 1;
            Im[0] = 0;
        }
        #endregion
    }
}
