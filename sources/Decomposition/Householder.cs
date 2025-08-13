using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines Householder transformation.
    /// <remarks>
    /// This is a linear transformation H (u) of the vector space V, which describes its mapping with respect to the hyperplane,
    /// which passes through the origin. It was proposed in 1958 by the American mathematician Elston Scott Householder. Widely used in linear algebra for QR decomposition of a matrix.
    /// In addition, the Householder transform is actively used for orthogonalization of bases; ultimately, the Householder matrix has the following properties:
    /// H = H', H' * H = I; det(H) = -1.
    /// In this class, two types of the Householder transform are implemented: reduction to a three-diagonal matrix and construction of the Householder matrix from a given vector.
    /// In the first case, the original Square matrix is defined as: A = H * T * H '.
    /// More information can be found on the website: 
    /// https://en.wikipedia.org/wiki/Householder_transformation
    /// </remarks>
    /// </summary>
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
        /// <param name="A">Square matrix</param>
        public Householder(float[,] A)
        {
            if (!Matrice.IsSymmetric(A))
                throw new Exception("The matrix must be symmetric");

            // properties:
            this.n = A.GetLength(0);

            Tridiagonal tridiagonal = new Tridiagonal(A);
            this.Re = tridiagonal.Diagonal;
            this.Im = tridiagonal.OffDiagonal;
            this.matrices = tridiagonal.Orthogonal;
        }
        #endregion

        #region Standart voids
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
        /// 
        /// </summary>
        /// <param name="v"></param>
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
        #endregion
    }
}
