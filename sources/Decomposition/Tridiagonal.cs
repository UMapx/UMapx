using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines tridiagonal decomposition of a symmetric matrix using Householder transformations.
    /// </summary>
    [Serializable]
    public class Tridiagonal
    {
        #region Private data
        private readonly int n;
        private readonly float[] d;
        private readonly float[] e;
        private readonly float[][] q;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes tridiagonal decomposition.
        /// </summary>
        /// <param name="A">Symmetric matrix</param>
        public Tridiagonal(float[,] A)
        {
            if (!Matrice.IsSymmetric(A))
                throw new Exception("The matrix must be symmetric");

            this.n = A.GetLength(0);
            this.d = new float[n];
            this.e = new float[n];
            this.q = Jagged.ToJagged(A);

            tred2();
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Gets orthogonal matrix.
        /// </summary>
        public float[,] Q
        {
            get { return Jagged.FromJagged(q); }
        }
        /// <summary>
        /// Gets tridiagonal matrix.
        /// </summary>
        public float[,] T
        {
            get
            {
                float[,] T = new float[n, n];
                int i;

                for (i = 0; i < n; i++)
                    T[i, i] = d[i];

                for (i = 1; i < n; i++)
                {
                    T[i - 1, i] = e[i];
                    T[i, i - 1] = e[i];
                }

                return T;
            }
        }
        #endregion

        #region Internal data
        internal float[] Diagonal => d;
        internal float[] OffDiagonal => e;
        internal float[][] Orthogonal => q;
        #endregion

        #region Private voids
        /// <summary>
        /// Symmetric Householder reduction to tridiagonal form.
        /// This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch, and Wilkinson,
        /// Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
        /// </summary>
        private void tred2()
        {
            int i, j, k;
            for (j = 0; j < n; j++)
                d[j] = q[n - 1][j];

            float scale, h, f, g, hh;

            // Householder reduction to tridiagonal form.
            for (i = n - 1; i > 0; i--)
            {
                // Scale to avoid under/overflow.
                scale = 0;
                h = 0;
                for (k = 0; k < i; k++)
                    scale += Math.Abs(d[k]);

                if (scale == 0)
                {
                    e[i] = d[i - 1];
                    for (j = 0; j < i; j++)
                    {
                        d[j] = q[i - 1][j];
                        q[i][j] = 0;
                        q[j][i] = 0;
                    }
                }
                else
                {
                    // Generate Householder vector.
                    for (k = 0; k < i; k++)
                    {
                        d[k] /= scale;
                        h += d[k] * d[k];
                    }

                    f = d[i - 1];
                    g = (float)Math.Sqrt(h);
                    if (f > 0) g = -g;

                    e[i] = scale * g;
                    h -= f * g;
                    d[i - 1] = f - g;
                    for (j = 0; j < i; j++)
                        e[j] = 0;

                    // Apply similarity transformation to remaining columns.
                    for (j = 0; j < i; j++)
                    {
                        f = d[j];
                        q[j][i] = f;
                        g = e[j] + q[j][j] * f;
                        for (k = j + 1; k <= i - 1; k++)
                        {
                            g += q[k][j] * d[k];
                            e[k] += q[k][j] * f;
                        }
                        e[j] = g;
                    }

                    f = 0;
                    for (j = 0; j < i; j++)
                    {
                        e[j] /= h;
                        f += e[j] * d[j];
                    }

                    hh = f / (h + h);
                    for (j = 0; j < i; j++)
                        e[j] -= hh * d[j];

                    for (j = 0; j < i; j++)
                    {
                        f = d[j];
                        g = e[j];
                        for (k = j; k <= i - 1; k++)
                            q[k][j] -= (f * e[k] + g * d[k]);

                        d[j] = q[i - 1][j];
                        q[i][j] = 0;
                    }
                }
                d[i] = h;
            }

            // Accumulate transformations.
            for (i = 0; i < n - 1; i++)
            {
                q[n - 1][i] = q[i][i];
                q[i][i] = 1;
                h = d[i + 1];
                if (h != 0)
                {
                    for (k = 0; k <= i; k++)
                        d[k] = q[k][i + 1] / h;

                    for (j = 0; j <= i; j++)
                    {
                        g = 0;
                        for (k = 0; k <= i; k++)
                            g += q[k][i + 1] * q[k][j];
                        for (k = 0; k <= i; k++)
                            q[k][j] -= g * d[k];
                    }
                }

                for (k = 0; k <= i; k++)
                    q[k][i + 1] = 0;
            }

            for (j = 0; j < n; j++)
            {
                d[j] = q[n - 1][j];
                q[n - 1][j] = 0;
            }

            q[n - 1][n - 1] = 1;
            e[0] = 0;
        }
        #endregion
    }
}
