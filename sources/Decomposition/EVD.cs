using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines eigenvalue decomposition.
    /// <remarks>
    /// The eigenvalue decomposition is the representation of the square matrix A in the form of the product of three matrices A = V * D * inv (V), where V is the matrix of spectral vectors and D is the diagonal (generally complex) matrix of eigenvalues.
    /// The matrix A can also be represented as the product of three matrices: A = V * R * inv (V), where R is a real almost diagonal eigenvalue matrix.
    /// Not all matrices can be represented in this form, but only those that have a complete set of eigenvectors.
    /// Eigenvalue decomposition can be used to find the eigenvalues ​​and eigenvectors of the matrix, solve linear systems of equations, invert the matrix, find the determinant of the matrix, and calculate the analytic functions of the matrices.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix
    /// </remarks>
    /// </summary>
    [Serializable]
    public class EVD
    {
        #region Private data
        private int n;
        private float[] Re, Im;
        private float[][] matrices;
        private float[][] hessenberg;
        private float[] orthogonal;
        private float eps;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes eigenvalue decomposition.
        /// </summary>
        /// <param name="A">Square matrix</param>
        /// <param name="eps">Epsilon [0, 1]</param>
        public EVD(float[,] A, float eps = 1e-16f)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("The matrix must be square");

            this.n = A.GetLength(0);
            this.Re = new float[n];
            this.Im = new float[n];
            this.eps = Maths.Float(eps);

            // for symmetric matrices eigen-value decomposition
            // without Hessenberg form.
            if (Matrice.IsSymmetric(A))
            {
                hessenberg = Jagged.Zero(n, n);
                matrices = Jagged.ToJagged(A);

                tred2(); // Tridiagonalize.
                tql2();  // Diagonalize.
            }
            // with Hessenberg form.
            else
            {
                matrices = Jagged.Zero(n, n);
                hessenberg = Jagged.ToJagged(A);
                orthogonal = new float[n];

                orthes(); // Reduce to Hessenberg form.
                hqr2();   // Reduce Hessenberg to real Schur form.
            }
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets eigenvectors.
        /// </summary>
        public float[,] V
        {
            get { return Jagged.FromJagged(matrices); }
        }
        /// <summary>
        /// Gets eigenvalues.
        /// </summary>
        public Complex32[] D
        {
            get
            {
                Complex32[] D = new Complex32[n];

                for (int i = 0; i < n; i++)
                {
                    D[i] = new Complex32(Re[i], Im[i]);
                }

                return D;
            }
        }
        /// <summary>
        /// Gets the real diagonal eigenvalue matrix.
        /// </summary>
        public float[,] R
        {
            get
            {
                float[,] D = new float[n, n];
                int i, j;

                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        D[i, j] = 0;
                    }

                    D[i, i] = Re[i];

                    if (Im[i] > 0)
                    {
                        D[i, i + 1] = Im[i];
                    }
                    else if (Im[i] < 0)
                    {
                        D[i, i - 1] = Im[i];
                    }
                }

                return D;
            }
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
        /// Symmetric Householder reduction to tridiagonal form.
        /// This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch, and Wilkinson, 
        /// Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK. 
        /// </summary>
        private void tred2()
        {
            int i, j, k;

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
                    scale = scale + Math.Abs(Re[k]);

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
        /// <summary>
        /// Symmetric tridiagonal QL algorithm.
        /// This is derived from the Algol procedures tql2, by Bowdler, Martin, Reinsch, and Wilkinson, 
        /// Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
        /// </summary>
        private void tql2()
        {
            float f = 0;
            float tst1 = 0;
            int i, l, j, k, iter, m;
            float g, p, r, dl1, h;
            float c, c2, c3, el1, s, s2;

            for (i = 1; i < n; i++)
                Im[i - 1] = Im[i];

            Im[n - 1] = 0;

            for (l = 0; l < n; l++)
            {
                // Find small subdiagonal element.
                tst1 = System.Math.Max(tst1, System.Math.Abs(Re[l]) + System.Math.Abs(Im[l]));
                m = l;
                while (m < n)
                {
                    if (System.Math.Abs(Im[m]) <= eps * tst1)
                        break;
                    m++;
                }

                // If m == l, d[l] is an eigenvalue, otherwise, iterate.
                if (m > l)
                {
                    iter = 0;
                    do
                    {
                        iter = iter + 1;  // (Could check iteration count here.)

                        // Compute implicit shift
                        g = Re[l];
                        p = (Re[l + 1] - g) / (2 * Im[l]);
                        r = Maths.Hypotenuse(p, 1);
                        if (p < 0)
                        {
                            r = -r;
                        }

                        Re[l] = Im[l] / (p + r);
                        Re[l + 1] = Im[l] * (p + r);
                        dl1 = Re[l + 1];
                        h = g - Re[l];
                        for (i = l + 2; i < n; i++)
                        {
                            Re[i] -= h;
                        }

                        f = f + h;

                        // Implicit QL transformation.
                        p = Re[m];
                        c = 1;
                        c2 = c;
                        c3 = c;
                        el1 = Im[l + 1];
                        s = 0;
                        s2 = 0;

                        for (i = m - 1; i >= l; i--)
                        {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * Im[i];
                            h = c * p;
                            r = Maths.Hypotenuse(p, Im[i]);
                            Im[i + 1] = s * r;
                            s = Im[i] / r;
                            c = p / r;
                            p = c * Re[i] - s * g;
                            Re[i + 1] = h + s * (c * g + s * Re[i]);

                            // Accumulate transformation.
                            for (k = 0; k < n; k++)
                            {
                                h = matrices[k][i + 1];
                                matrices[k][i + 1] = s * matrices[k][i] + c * h;
                                matrices[k][i] = c * matrices[k][i] - s * h;
                            }
                        }

                        p = -s * s2 * c3 * el1 * Im[l] / dl1;
                        Im[l] = s * p;
                        Re[l] = c * p;

                        // Check for convergence.
                    }
                    while (System.Math.Abs(Im[l]) > eps * tst1);
                }
                Re[l] = Re[l] + f;
                Im[l] = 0;
            }

            // Sort eigenvalues and corresponding Matrices.
            for (i = 0; i < n - 1; i++)
            {
                k = i;
                p = Re[i];
                for (j = i + 1; j < n; j++)
                {
                    if (Re[j] < p)
                    {
                        k = j;
                        p = Re[j];
                    }
                }

                if (k != i)
                {
                    Re[k] = Re[i];
                    Re[i] = p;
                    for (j = 0; j < n; j++)
                    {
                        p = matrices[j][i];
                        matrices[j][i] = matrices[j][k];
                        matrices[j][k] = p;
                    }
                }
            }
        }
        /// <summary>
        /// Nonsymmetric reduction to Hessenberg form.
        /// This is derived from the Algol procedures orthes and ortran, by Martin and Wilkinson, 
        /// Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutines in EISPACK.
        /// </summary>
        private void orthes()
        {
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
        }
        /// <summary>
        /// Nonsymmetric reduction from Hessenberg to real Schur form.   
        /// This is derived from the Algol procedure hqr2, by Martin and Wilkinson, Handbook for Auto. Comp.,
        /// Vol.ii-Linear Algebra, and the corresponding  Fortran subroutine in EISPACK.
        /// </summary>
        private void hqr2()
        {
            int nn = this.n;
            int n = nn - 1;
            int low = 0;
            int high = nn - 1;
            //float eps = 2 * float.Epsilon;
            float exshift = 0;
            float p = 0;
            float q = 0;
            float r = 0;
            float s = 0;
            float z = 0;
            float t;
            float w;
            float x;
            float y;
            int i, j, k, m;
            bool notlast;

            // Store roots isolated by balanc and compute matrix norm
            float norm = 0;
            for (i = 0; i < nn; i++)
            {
                if (i < low | i > high)
                {
                    Re[i] = hessenberg[i][i];
                    Im[i] = 0;
                }

                for (j = System.Math.Max(i - 1, 0); j < nn; j++)
                    norm = norm + System.Math.Abs(hessenberg[i][j]);
            }

            // Outer loop over eigenvalue index
            int iter = 0;
            while (n >= low)
            {
                // Look for single small sub-diagonal element
                int l = n;
                while (l > low)
                {
                    s = System.Math.Abs(hessenberg[l - 1][l - 1]) + System.Math.Abs(hessenberg[l][l]);

                    if (s == 0)
                        s = norm;

                    if (float.IsNaN(s))
                        break;

                    if (System.Math.Abs(hessenberg[l][l - 1]) < eps * s)
                        break;

                    l--;
                }

                // Check for convergence
                if (l == n)
                {
                    // One root found
                    hessenberg[n][n] = hessenberg[n][n] + exshift;
                    Re[n] = hessenberg[n][n];
                    Im[n] = 0;
                    n--;
                    iter = 0;
                }
                else if (l == n - 1)
                {
                    // Two roots found
                    w = hessenberg[n][n - 1] * hessenberg[n - 1][n];
                    p = (hessenberg[n - 1][n - 1] - hessenberg[n][n]) / 2;
                    q = p * p + w;
                    z = (float)System.Math.Sqrt(System.Math.Abs(q));
                    hessenberg[n][n] = hessenberg[n][n] + exshift;
                    hessenberg[n - 1][n - 1] = hessenberg[n - 1][n - 1] + exshift;
                    x = hessenberg[n][n];

                    if (q >= 0)
                    {
                        // Real pair
                        z = (p >= 0) ? (p + z) : (p - z);
                        Re[n - 1] = x + z;
                        Re[n] = Re[n - 1];
                        if (z != 0)
                            Re[n] = x - w / z;
                        Im[n - 1] = 0;
                        Im[n] = 0;
                        x = hessenberg[n][n - 1];
                        s = System.Math.Abs(x) + System.Math.Abs(z);
                        p = x / s;
                        q = z / s;
                        r = (float)System.Math.Sqrt(p * p + q * q);
                        p = p / r;
                        q = q / r;

                        // Row modification
                        for (j = n - 1; j < nn; j++)
                        {
                            z = hessenberg[n - 1][j];
                            hessenberg[n - 1][j] = q * z + p * hessenberg[n][j];
                            hessenberg[n][j] = q * hessenberg[n][j] - p * z;
                        }

                        // Column modification
                        for (i = 0; i <= n; i++)
                        {
                            z = hessenberg[i][n - 1];
                            hessenberg[i][n - 1] = q * z + p * hessenberg[i][n];
                            hessenberg[i][n] = q * hessenberg[i][n] - p * z;
                        }

                        // Accumulate transformations
                        for (i = low; i <= high; i++)
                        {
                            z = matrices[i][n - 1];
                            matrices[i][n - 1] = q * z + p * matrices[i][n];
                            matrices[i][n] = q * matrices[i][n] - p * z;
                        }
                    }
                    else
                    {
                        // Complex pair
                        Re[n - 1] = x + p;
                        Re[n] = x + p;
                        Im[n - 1] = z;
                        Im[n] = -z;
                    }

                    n = n - 2;
                    iter = 0;
                }
                else
                {
                    // No convergence yet     

                    // Form shift
                    x = hessenberg[n][n];
                    y = 0;
                    w = 0;
                    if (l < n)
                    {
                        y = hessenberg[n - 1][n - 1];
                        w = hessenberg[n][n - 1] * hessenberg[n - 1][n];
                    }

                    // Wilkinson's original ad hoc shift
                    if (iter == 10)
                    {
                        exshift += x;
                        for (i = low; i <= n; i++)
                            hessenberg[i][i] -= x;

                        s = System.Math.Abs(hessenberg[n][n - 1]) + System.Math.Abs(hessenberg[n - 1][n - 2]);
                        x = y = (float)0.75 * s;
                        w = (float)(-0.4375) * s * s;
                    }

                    // MATLAB's new ad hoc shift
                    if (iter == 30)
                    {
                        s = (y - x) / 2;
                        s = s * s + w;
                        if (s > 0)
                        {
                            s = (float)System.Math.Sqrt(s);
                            if (y < x) s = -s;
                            s = x - w / ((y - x) / 2 + s);
                            for (i = low; i <= n; i++)
                                hessenberg[i][i] -= s;
                            exshift += s;
                            x = y = w = (float)0.964;
                        }
                    }

                    iter = iter + 1;

                    // Look for two consecutive small sub-diagonal elements
                    m = n - 2;
                    while (m >= l)
                    {
                        z = hessenberg[m][m];
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / hessenberg[m + 1][m] + hessenberg[m][m + 1];
                        q = hessenberg[m + 1][m + 1] - z - r - s;
                        r = hessenberg[m + 2][m + 1];
                        s = System.Math.Abs(p) + System.Math.Abs(q) + System.Math.Abs(r);
                        p = p / s;
                        q = q / s;
                        r = r / s;
                        if (m == l)
                            break;
                        if (System.Math.Abs(hessenberg[m][m - 1]) * (System.Math.Abs(q) + System.Math.Abs(r)) < eps * (System.Math.Abs(p) * (System.Math.Abs(hessenberg[m - 1][m - 1]) + System.Math.Abs(z) + System.Math.Abs(hessenberg[m + 1][m + 1]))))
                            break;
                        m--;
                    }

                    for (i = m + 2; i <= n; i++)
                    {
                        hessenberg[i][i - 2] = 0;
                        if (i > m + 2)
                            hessenberg[i][i - 3] = 0;
                    }

                    // float QR step involving rows l:n and columns m:n
                    for (k = m; k <= n - 1; k++)
                    {
                        notlast = (k != n - 1);
                        if (k != m)
                        {
                            p = hessenberg[k][k - 1];
                            q = hessenberg[k + 1][k - 1];
                            r = (notlast ? hessenberg[k + 2][k - 1] : 0);
                            x = System.Math.Abs(p) + System.Math.Abs(q) + System.Math.Abs(r);
                            if (x != 0)
                            {
                                p = p / x;
                                q = q / x;
                                r = r / x;
                            }
                        }

                        if (x == 0) break;

                        s = (float)System.Math.Sqrt(p * p + q * q + r * r);
                        if (p < 0) s = -s;

                        if (s != 0)
                        {
                            if (k != m)
                                hessenberg[k][k - 1] = -s * x;
                            else
                                if (l != m)
                                hessenberg[k][k - 1] = -hessenberg[k][k - 1];

                            p = p + s;
                            x = p / s;
                            y = q / s;
                            z = r / s;
                            q = q / p;
                            r = r / p;

                            // Row modification
                            for (j = k; j < nn; j++)
                            {
                                p = hessenberg[k][j] + q * hessenberg[k + 1][j];
                                if (notlast)
                                {
                                    p = p + r * hessenberg[k + 2][j];
                                    hessenberg[k + 2][j] = hessenberg[k + 2][j] - p * z;
                                }

                                hessenberg[k][j] = hessenberg[k][j] - p * x;
                                hessenberg[k + 1][j] = hessenberg[k + 1][j] - p * y;
                            }

                            // Column modification
                            for (i = 0; i <= System.Math.Min(n, k + 3); i++)
                            {
                                p = x * hessenberg[i][k] + y * hessenberg[i][k + 1];
                                if (notlast)
                                {
                                    p = p + z * hessenberg[i][k + 2];
                                    hessenberg[i][k + 2] = hessenberg[i][k + 2] - p * r;
                                }

                                hessenberg[i][k] = hessenberg[i][k] - p;
                                hessenberg[i][k + 1] = hessenberg[i][k + 1] - p * q;
                            }

                            // Accumulate transformations
                            for (i = low; i <= high; i++)
                            {
                                p = x * matrices[i][k] + y * matrices[i][k + 1];
                                if (notlast)
                                {
                                    p = p + z * matrices[i][k + 2];
                                    matrices[i][k + 2] = matrices[i][k + 2] - p * r;
                                }

                                matrices[i][k] = matrices[i][k] - p;
                                matrices[i][k + 1] = matrices[i][k + 1] - p * q;
                            }
                        }
                    }
                }
            }

            // Backsubstitute to find Matrices of upper triangular form
            if (norm == 0)
            {
                return;
            }

            for (n = nn - 1; n >= 0; n--)
            {
                p = Re[n];
                q = Im[n];

                // Real Matrice
                if (q == 0)
                {
                    int l = n;
                    hessenberg[n][n] = 1;
                    for (i = n - 1; i >= 0; i--)
                    {
                        w = hessenberg[i][i] - p;
                        r = 0;
                        for (j = l; j <= n; j++)
                            r = r + hessenberg[i][j] * hessenberg[j][n];

                        if (Im[i] < 0)
                        {
                            z = w;
                            s = r;
                        }
                        else
                        {
                            l = i;
                            if (Im[i] == 0)
                            {
                                hessenberg[i][n] = (w != 0) ? (-r / w) : (-r / (eps * norm));
                            }
                            else
                            {
                                // Solve real equations
                                x = hessenberg[i][i + 1];
                                y = hessenberg[i + 1][i];
                                q = (Re[i] - p) * (Re[i] - p) + Im[i] * Im[i];
                                t = (x * s - z * r) / q;
                                hessenberg[i][n] = t;
                                hessenberg[i + 1][n] = (System.Math.Abs(x) > System.Math.Abs(z)) ? ((-r - w * t) / x) : ((-s - y * t) / z);
                            }

                            // Overflow control
                            t = System.Math.Abs(hessenberg[i][n]);
                            if ((eps * t) * t > 1)
                                for (j = i; j <= n; j++)
                                    hessenberg[j][n] = hessenberg[j][n] / t;
                        }
                    }
                }
                else if (q < 0)
                {
                    // Complex Matrice
                    int l = n - 1;

                    // Last Matrice component imaginary so matrix is triangular
                    if (System.Math.Abs(hessenberg[n][n - 1]) > System.Math.Abs(hessenberg[n - 1][n]))
                    {
                        hessenberg[n - 1][n - 1] = q / hessenberg[n][n - 1];
                        hessenberg[n - 1][n] = -(hessenberg[n][n] - p) / hessenberg[n][n - 1];
                    }
                    else
                    {
                        cdiv(0, -hessenberg[n - 1][n], hessenberg[n - 1][n - 1] - p, q, ref hessenberg[n - 1][n - 1], ref hessenberg[n - 1][n]);
                    }

                    hessenberg[n][n - 1] = 0;
                    hessenberg[n][n] = 1;
                    for (i = n - 2; i >= 0; i--)
                    {
                        float ra, sa, vr, vi;
                        ra = 0;
                        sa = 0;
                        for (j = l; j <= n; j++)
                        {
                            ra = ra + hessenberg[i][j] * hessenberg[j][n - 1];
                            sa = sa + hessenberg[i][j] * hessenberg[j][n];
                        }

                        w = hessenberg[i][i] - p;

                        if (Im[i] < 0)
                        {
                            z = w;
                            r = ra;
                            s = sa;
                        }
                        else
                        {
                            l = i;
                            if (Im[i] == 0)
                            {
                                cdiv(-ra, -sa, w, q, ref hessenberg[i][n - 1], ref hessenberg[i][n]);
                            }
                            else
                            {
                                // Solve complex equations
                                x = hessenberg[i][i + 1];
                                y = hessenberg[i + 1][i];
                                vr = (Re[i] - p) * (Re[i] - p) + Im[i] * Im[i] - q * q;
                                vi = (Re[i] - p) * 2 * q;
                                if (vr == 0 & vi == 0)
                                    vr = eps * norm * (System.Math.Abs(w) + System.Math.Abs(q) + System.Math.Abs(x) + System.Math.Abs(y) + System.Math.Abs(z));
                                cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi, ref hessenberg[i][n - 1], ref hessenberg[i][n]);
                                if (System.Math.Abs(x) > (System.Math.Abs(z) + System.Math.Abs(q)))
                                {
                                    hessenberg[i + 1][n - 1] = (-ra - w * hessenberg[i][n - 1] + q * hessenberg[i][n]) / x;
                                    hessenberg[i + 1][n] = (-sa - w * hessenberg[i][n] - q * hessenberg[i][n - 1]) / x;
                                }
                                else
                                {
                                    cdiv(-r - y * hessenberg[i][n - 1], -s - y * hessenberg[i][n], z, q, ref hessenberg[i + 1][n - 1], ref hessenberg[i + 1][n]);
                                }
                            }

                            // Overflow control
                            t = System.Math.Max(System.Math.Abs(hessenberg[i][n - 1]), System.Math.Abs(hessenberg[i][n]));
                            if ((eps * t) * t > 1)
                            {
                                for (j = i; j <= n; j++)
                                {
                                    hessenberg[j][n - 1] = hessenberg[j][n - 1] / t;
                                    hessenberg[j][n] = hessenberg[j][n] / t;
                                }
                            }
                        }
                    }
                }
            }

            // Matrices of isolated roots
            for (i = 0; i < nn; i++)
                if (i < low | i > high)
                    for (j = i; j < nn; j++)
                        matrices[i][j] = hessenberg[i][j];

            // Back transformation to get eigenMatrices of original matrix
            for (j = nn - 1; j >= low; j--)
            {
                for (i = low; i <= high; i++)
                {
                    z = 0;
                    for (k = low; k <= System.Math.Min(j, high); k++)
                        z = z + matrices[i][k] * hessenberg[k][j];
                    matrices[i][j] = z;
                }
            }
        }
        /// <summary>
        /// Complex scalar division using a numerically stable branch (Smith’s method).
        /// Computes (xr + i·xi) / (yr + i·yi) and stores the real/imag parts in <paramref name="cdivr"/> / <paramref name="cdivi"/>.
        /// </summary>
        /// <param name="xr">Real part of the numerator.</param>
        /// <param name="xi">Imag part of the numerator.</param>
        /// <param name="yr">Real part of the denominator.</param>
        /// <param name="yi">Imag part of the denominator.</param>
        /// <param name="cdivr">[out] Real part of the quotient.</param>
        /// <param name="cdivi">[out] Imag part of the quotient.</param>
        /// <remarks>
        /// Chooses the scaling branch by comparing |yr| and |yi| to avoid overflow/underflow.
        /// If both <paramref name="yr"/> and <paramref name="yi"/> are zero, the result follows IEEE-754 (Inf/NaN).
        /// </remarks>
        private static void cdiv(float xr, float xi, float yr, float yi, ref float cdivr, ref float cdivi)
        {
            // Complex scalar division.
            float r;
            float d;

            if (System.Math.Abs(yr) > System.Math.Abs(yi))
            {
                r = yi / yr;
                d = yr + r * yi;
                cdivr = (xr + r * xi) / d;
                cdivi = (xi - r * xr) / d;
            }
            else
            {
                r = yr / yi;
                d = yi + r * yr;
                cdivr = (r * xr + xi) / d;
                cdivi = (r * xi - xr) / d;
            }
        }

        #endregion
    }
}
