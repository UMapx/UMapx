using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines Schur decomposition.
    /// </summary>
    /// <remarks>
    /// This is a representation of a square matrix in the form of a product of three matrices: A = Q * T * Qᵀ,
    /// where Q is a unitary matrix and T is a quasi upper triangular matrix (Schur form).
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Schur_decomposition
    /// </remarks>
    [Serializable]
    public class Schur
    {
        #region Private data
        private int n;
        private float[][] matrices; // unitary matrix
        private float[][] hessenberg; // schur form
        private double[] Re, Im;
        private float eps;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes Schur decomposition.
        /// </summary>
        /// <param name="A">Nonempty square matrix with finite real entries.</param>
        /// <param name="eps">Relative deflation tolerance, clamped to [0, 1] with a double-roundoff floor.</param>
        /// <exception cref="ArgumentException">The input is empty, nonsquare, or contains nonfinite entries.</exception>
        /// <exception cref="InvalidOperationException">The QR iteration limit is reached before convergence.</exception>
        public Schur(float[,] A, float eps = 1e-16f)
        {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (A.GetLength(0) == 0) throw new ArgumentException("The matrix must be nonempty.", nameof(A));
            if (float.IsNaN(eps)) throw new ArgumentOutOfRangeException(nameof(eps));
            if (!Matrice.IsSquare(A))
                throw new ArgumentException("The matrix must be square");

            this.n = A.GetLength(0);
            this.Re = new double[n];
            this.Im = new double[n];
            this.eps = Maths.Float(eps);

            var hessenberg = ScaleInput(A, out double inputScale);
            var matrices = ReduceToHessenberg(hessenberg);
            this.matrices = Jagged.Zero(n, n);
            this.hessenberg = Jagged.Zero(n, n);
            hqr2(hessenberg, matrices, inputScale);
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the unitary matrix Q.
        /// </summary>
        public float[,] Q
        {
            get { return Jagged.FromJagged(matrices); }
        }
        /// <summary>
        /// Gets the quasi upper triangular matrix T.
        /// </summary>
        public float[,] T
        {
            get { return Jagged.FromJagged(hessenberg); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Scales a finite square matrix before Hessenberg reduction to protect products in the QR shifts.
        /// </summary>
        /// <param name="matrix">Original real square matrix.</param>
        /// <param name="scale">Receives the maximum magnitude, or one for a zero matrix.</param>
        /// <returns>A scaled copy; the Schur form must subsequently be multiplied by scale.</returns>
        private static double[][] ScaleInput(float[,] matrix, out double scale)
        {
            scale = 0;
            foreach (float value in matrix)
            {
                if (float.IsNaN(value) || float.IsInfinity(value))
                    throw new ArgumentException("The matrix must contain only finite values.", nameof(matrix));
                scale = Math.Max(scale, Math.Abs((double)value));
            }
            if (scale == 0) scale = 1;
            var result = new double[matrix.GetLength(0)][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[matrix.GetLength(1)];
                for (int j = 0; j < result[i].Length; j++) result[i][j] = matrix[i, j] / scale;
            }
            return result;
        }
        /// <summary>
        /// Reduces a scaled square matrix to Hessenberg form by Householder similarities in double precision.
        /// </summary>
        /// <param name="hessenberg">Matrix overwritten by its upper Hessenberg form H.</param>
        /// <returns>The orthogonal accumulator Q satisfying A = Q*H*Q^T before QR iteration.</returns>
        private static double[][] ReduceToHessenberg(double[][] hessenberg)
        {
            int n = hessenberg.Length;
            var matrices = new double[n][];
            for (int row = 0; row < n; row++) matrices[row] = new double[n];
            var orthogonal = new double[n];
            int low = 0;
            int high = n - 1;
            int m, i, j;
            double scale, h, g, f;

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

                    g = System.Math.Sqrt(h);
                    if (orthogonal[m] > 0) g = -g;

                    h = h - orthogonal[m] * g;
                    orthogonal[m] = orthogonal[m] - g;

                    // Apply Householder similarity transformation
                    // H = (I - u * u' / h) * H * (I - u * u' / h).
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

                        // Double division avoids possible underflow.
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

            return matrices;
        }
        /// <summary>
        /// Reduces Hessenberg form to real Schur form using bounded double-shift QR iteration.
        /// </summary>
        /// <param name="hessenberg">Scaled Hessenberg matrix, overwritten by its Schur form.</param>
        /// <param name="matrices">Orthogonal reduction accumulator, overwritten by the Schur vectors.</param>
        /// <param name="inputScale">Positive input scale restored when storing the single-precision Schur form.</param>
        private void hqr2(double[][] hessenberg, double[][] matrices, double inputScale)
        {
            int nn = this.n;
            double eps = Math.Max(this.eps, 2.2204460492503131e-16);
            int n = nn - 1;
            int low = 0;
            int high = nn - 1;
            double exshift = 0;
            double p = 0;
            double q = 0;
            double r = 0;
            double s = 0;
            double z = 0;
            double w;
            double x;
            double y;
            int i, j, k, m;
            bool notlast;

            // Store roots isolated by balanc and compute matrix norm
            double norm = 0;
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
                    // Exact zeros must deflate even when the matrix norm or eps is zero.
                    if (System.Math.Abs(hessenberg[l][l - 1]) <= eps * s)
                    {
                        hessenberg[l][l - 1] = 0;
                        break;
                    }
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
                    z = System.Math.Sqrt(System.Math.Abs(q));
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
                        r = System.Math.Sqrt(p * p + q * q);
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
                        for (i = 0; i < nn; i++)
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
                    if (Im[n] == 0) hessenberg[n][n - 1] = 0;
                    n = n - 2;
                    iter = 0;
                }
                else
                {
                    // No convergence yet
                    x = hessenberg[n][n];
                    y = 0;
                    w = 0;
                    if (l < n)
                    {
                        y = hessenberg[n - 1][n - 1];
                        w = hessenberg[n][n - 1] * hessenberg[n - 1][n];
                    }

                    if (iter == 10)
                    {
                        exshift += x;
                        for (i = low; i <= n; i++)
                            hessenberg[i][i] -= x;
                        s = System.Math.Abs(hessenberg[n][n - 1]) + System.Math.Abs(hessenberg[n - 1][n - 2]);
                        x = y = 0.75f * s;
                        w = -0.4375f * s * s;
                    }

                    if (iter == 30)
                    {
                        s = (y - x) / 2;
                        s = s * s + w;
                        if (s > 0)
                        {
                            s = System.Math.Sqrt(s);
                            if (y < x)
                                s = -s;
                            s = x - w / ((y - x) / 2 + s);
                            for (i = low; i <= n; i++)
                                hessenberg[i][i] -= s;
                            exshift += s;
                            x = y = w = 0.964f;
                        }
                    }

                    // A failed iteration must not leave the caller in an unbounded loop.
                    if (++iter > 100 * nn)
                        throw new InvalidOperationException("Schur decomposition failed to converge.");

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
                        if (System.Math.Abs(hessenberg[m][m - 1]) * (System.Math.Abs(q) + System.Math.Abs(r)) <
                            eps * (System.Math.Abs(p) * (System.Math.Abs(hessenberg[m - 1][m - 1]) + System.Math.Abs(z) + System.Math.Abs(hessenberg[m + 1][m + 1]))))
                            break;
                        m--;
                    }

                    for (i = m + 2; i <= n; i++)
                    {
                        hessenberg[i][i - 2] = 0;
                        if (i > m + 2)
                            hessenberg[i][i - 3] = 0;
                    }

                    // Double-shift QR step involving rows l:n and columns m:n
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

                        if (k != m && x == 0)
                            break;

                        s = System.Math.Sqrt(p * p + q * q + r * r);
                        if (p < 0)
                            s = -s;

                        if (s != 0)
                        {
                            if (k != m)
                                hessenberg[k][k - 1] = -s * x;
                            else if (l != m)
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
                            for (i = 0; i < nn; i++)
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
            for (int row = 0; row < nn; row++)
                for (int column = 0; column < nn; column++)
                {
                    // Entries below the first subdiagonal are structural zeros;
                    // roundoff from rotations must not leak into the public Schur form.
                    this.hessenberg[row][column] = column < row - 1 ? 0 : (float)(hessenberg[row][column] * inputScale);
                    this.matrices[row][column] = (float)matrices[row][column];
                }
        }
        #endregion
    }
}
