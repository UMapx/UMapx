using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines singular value decomposition.
    /// </summary>
    /// <remarks>
    /// This is a representation of a rectangular matrix A in the form of the product of three matrices A = U * S * Vᵀ, 
    /// where U are left vectors, V are right vectors, and S are singular values.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Singular_value_decomposition
    /// </remarks>
    [Serializable]
    public class SVD
    {
        #region Private data
        private int n, m;
        private int iterations;
        private float[][] Ur;
        private float[][] Vr;
        private float[] Sr;
        private bool reversed;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes singular value decomposition.
        /// </summary>
        /// <param name="A">Nonempty rectangular matrix with finite real entries.</param>
        /// <param name="iterations">Positive maximum number of QR sweeps per singular value.</param>
        /// <exception cref="ArgumentException">The matrix is empty or contains nonfinite entries.</exception>
        /// <exception cref="InvalidOperationException">The QR iteration limit is reached before convergence.</exception>
        public SVD(float[,] A, int iterations = 10)
        {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (A.GetLength(0) == 0 || A.GetLength(1) == 0)
                throw new ArgumentException("The matrix must be nonempty.", nameof(A));
            if (iterations < 1) throw new ArgumentOutOfRangeException(nameof(iterations), "The iteration limit must be positive.");
            // set:
            this.iterations = iterations;
            this.n = A.GetLength(0);
            this.m = A.GetLength(1);

            // options:
            if (n < m)
            {
                this.reversed = true;
                this.n = A.GetLength(1);
                this.m = A.GetLength(0);
                this.svdcmp(A.Transpose());
            }
            else
            {
                this.reversed = false;
                this.svdcmp(A);
            }
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the left vectors.
        /// </summary>
        public float[,] U
        {
            get
            {
                return reversed ? Jagged.FromJagged(Vr) : Jagged.FromJagged(Ur);
            }
        }
        /// <summary>
        /// Gets singular values.
        /// </summary>
        public float[] S
        {
            get { return Sr; }
        }
        /// <summary>
        /// Gets the right vectors.
        /// </summary>
        public float[,] V
        {
            get
            {
                return reversed ? Jagged.FromJagged(Ur) : Jagged.FromJagged(Vr);
            }
        }
        /// <summary>
        /// Gets the pseudoinverse matrix.
        /// </summary>
        /// <remarks>
        /// Singular values at or below max(rows, columns) * 2^-23 * max(S) are
        /// treated as numerically zero. This rank threshold uses single-precision
        /// input resolution; retained singular values are inverted in double precision.
        /// </remarks>
        public float[,] P
        {
            get
            {
                int rows = reversed ? m : n, columns = reversed ? n : m;
                var result = new float[columns, rows];
                double largest = 0;
                foreach (float singular in Sr) largest = Math.Max(largest, singular);
                double threshold = Math.Max(n, m) * 1.1920928955078125e-7 * largest;
                for (int i = 0; i < columns; i++)
                    for (int j = 0; j < rows; j++)
                    {
                        double sum = 0;
                        for (int k = 0; k < m; k++)
                            if (Sr[k] > threshold)
                                // Divide before narrowing to float, so a small singular value
                                // cannot create infinity times a zero vector component.
                                sum += reversed ? (double)Ur[i][k] * Vr[j][k] / Sr[k]
                                                : (double)Vr[i][k] * Ur[j][k] / Sr[k];
                        result[i, j] = (float)sum;
                    }
                return result;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Core SVD routine with double-precision work buffers for real single-precision inputs.
        /// Performs Householder bidiagonalization followed by Golub–Kahan QR iterations
        /// to compute singular values and left/right singular vectors.
        /// Populates the private fields: <c>Ur</c> (left vectors), <c>Vr</c> (right vectors),
        /// and <c>Sr</c> (non-negative singular values).
        /// </summary>
        /// <param name="A">
        /// Input matrix of size n×m. Assumes n ≥ m when called (the caller transposes
        /// beforehand if needed). The method works on an internal copy (jagged buffers)
        /// </param>
        /// <remarks>
        /// Uses jagged arrays for speed. Columns of U and V are orthonormal. The number
        /// of QR sweeps is limited by the instance field <see cref="iterations"/>.
        /// </remarks>
        private void svdcmp(float[,] A)
        {
            var Ur = ScaleInput(A, out double inputScale);
            var Sr = new double[m];
            var Vr = new double[m][];
            for (int row = 0; row < m; row++) Vr[row] = new double[m];
            double[] rv1 = new double[m];

            int flag, i, its, j, jj, k, l = 0, nm = 0;
            double anorm, c, f, g, h, e, scale, x, y, z;

            // householder reduction to bidiagonal form
            g = scale = anorm = 0.0f;

            for (i = 0; i < m; i++)
            {
                l = i + 1;
                rv1[i] = scale * g;
                g = e = scale = 0;

                if (i < n)
                {
                    for (k = i; k < n; k++)
                    {
                        scale += Math.Abs(Ur[k][i]);
                    }

                    if (scale != 0.0)
                    {
                        for (k = i; k < n; k++)
                        {
                            Ur[k][i] /= scale;
                            e += Ur[k][i] * Ur[k][i];
                        }

                        f = Ur[i][i];
                        g = -CopySign(Math.Sqrt(e), f);
                        h = f * g - e;
                        Ur[i][i] = f - g;

                        if (i != m - 1)
                        {
                            for (j = l; j < m; j++)
                            {
                                for (e = 0.0f, k = i; k < n; k++)
                                {
                                    e += Ur[k][i] * Ur[k][j];
                                }

                                f = e / h;

                                for (k = i; k < n; k++)
                                {
                                    Ur[k][j] += f * Ur[k][i];
                                }
                            }
                        }

                        for (k = i; k < n; k++)
                        {
                            Ur[k][i] *= scale;
                        }
                    }
                }

                Sr[i] = scale * g;
                g = e = scale = 0.0f;

                if ((i < n) && (i != m - 1))
                {
                    for (k = l; k < m; k++)
                    {
                        scale += Math.Abs(Ur[i][k]);
                    }

                    if (scale != 0.0)
                    {
                        for (k = l; k < m; k++)
                        {
                            Ur[i][k] /= scale;
                            e += Ur[i][k] * Ur[i][k];
                        }

                        f = Ur[i][l];
                        g = -CopySign(Math.Sqrt(e), f);
                        h = f * g - e;
                        Ur[i][l] = f - g;

                        for (k = l; k < m; k++)
                        {
                            rv1[k] = Ur[i][k] / h;
                        }

                        if (i != n - 1)
                        {
                            for (j = l; j < n; j++)
                            {
                                for (e = 0.0f, k = l; k < m; k++)
                                {
                                    e += Ur[j][k] * Ur[i][k];
                                }
                                for (k = l; k < m; k++)
                                {
                                    Ur[j][k] += e * rv1[k];
                                }
                            }
                        }

                        for (k = l; k < m; k++)
                        {
                            Ur[i][k] *= scale;
                        }
                    }
                }
                anorm = Math.Max(anorm, (Math.Abs(Sr[i]) + Math.Abs(rv1[i])));
            }

            // accumulation of right-hand transformations
            for (i = m - 1; i >= 0; i--)
            {
                if (i < m - 1)
                {
                    if (g != 0.0)
                    {
                        for (j = l; j < m; j++)
                        {
                            Vr[j][i] = (Ur[i][j] / Ur[i][l]) / g;
                        }

                        for (j = l; j < m; j++)
                        {
                            for (e = 0, k = l; k < m; k++)
                            {
                                e += Ur[i][k] * Vr[k][j];
                            }
                            for (k = l; k < m; k++)
                            {
                                Vr[k][j] += e * Vr[k][i];
                            }
                        }
                    }
                    for (j = l; j < m; j++)
                    {
                        Vr[i][j] = Vr[j][i] = 0;
                    }
                }
                Vr[i][i] = 1;
                g = rv1[i];
                l = i;
            }

            // accumulation of left-hand transformations
            for (i = m - 1; i >= 0; i--)
            {
                l = i + 1;
                g = Sr[i];

                if (i < m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        Ur[i][j] = 0.0f;
                    }
                }

                if (g != 0)
                {
                    g = 1.0f / g;

                    if (i != m - 1)
                    {
                        for (j = l; j < m; j++)
                        {
                            for (e = 0, k = l; k < n; k++)
                            {
                                e += Ur[k][i] * Ur[k][j];
                            }

                            f = (e / Ur[i][i]) * g;

                            for (k = i; k < n; k++)
                            {
                                Ur[k][j] += f * Ur[k][i];
                            }
                        }
                    }

                    for (j = i; j < n; j++)
                    {
                        Ur[j][i] *= g;
                    }
                }
                else
                {
                    for (j = i; j < n; j++)
                    {
                        Ur[j][i] = 0;
                    }
                }
                ++Ur[i][i];
            }

            // diagonalization of the bidiagonal form: Loop over singular values
            // and over allowed iterations
            for (k = m - 1; k >= 0; k--)
            {
                for (its = 0; its <= iterations; its++)
                {
                    flag = 1;

                    for (l = k; l >= 0; l--)
                    {
                        // test for splitting
                        nm = l - 1;

                        if (Math.Abs(rv1[l]) + anorm == anorm)
                        {
                            flag = 0;
                            break;
                        }

                        if (Math.Abs(Sr[nm]) + anorm == anorm)
                            break;
                    }

                    if (flag != 0)
                    {
                        c = 0.0f;
                        e = 1.0f;
                        for (i = l; i <= k; i++)
                        {
                            f = e * rv1[i];
                            rv1[i] *= c;

                            if (Math.Abs(f) + anorm == anorm) break;
                            g = Sr[i];
                            h = Hypotenuse(f, g);
                            Sr[i] = h;
                            h = 1.0f / h;
                            c = g * h;
                            e = -f * h;

                            // Apply the cancellation rotation to every row, including row zero.
                            for (j = 0; j < n; j++)
                            {
                                y = Ur[j][nm];
                                z = Ur[j][i];
                                Ur[j][nm] = y * c + z * e;
                                Ur[j][i] = z * c - y * e;
                            }
                        }
                    }

                    z = Sr[k];

                    if (l == k)
                    {
                        // convergence
                        if (z < 0.0)
                        {
                            // singular value is made nonnegative
                            Sr[k] = -z;

                            for (j = 0; j < m; j++)
                            {
                                Vr[j][k] = -Vr[j][k];
                            }
                        }
                        break;
                    }

                    if (its == iterations)
                        throw new InvalidOperationException("Singular value decomposition failed to converge within the iteration limit.");

                    // shift from bottom 2-by-2 minor
                    x = Sr[l];
                    nm = k - 1;
                    y = Sr[nm];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);
                    g = Hypotenuse(f, 1.0f);
                    f = ((x - z) * (x + z) + h * ((y / (f + CopySign(g, f))) - h)) / x;

                    // next QR transformation
                    c = e = 1.0f;

                    for (j = l; j <= nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = Sr[i];
                        h = e * g;
                        g = c * g;
                        z = Hypotenuse(f, h);
                        rv1[j] = z;
                        c = f / z;
                        e = h / z;
                        f = x * c + g * e;
                        g = g * c - x * e;
                        h = y * e;
                        y *= c;

                        for (jj = 0; jj < m; jj++)
                        {
                            x = Vr[jj][j];
                            z = Vr[jj][i];
                            Vr[jj][j] = x * c + z * e;
                            Vr[jj][i] = z * c - x * e;
                        }

                        z = Hypotenuse(f, h);
                        Sr[j] = z;

                        if (z != 0)
                        {
                            z = 1.0f / z;
                            c = f * z;
                            e = h * z;
                        }

                        f = c * g + e * y;
                        x = c * y - e * g;

                        for (jj = 0; jj < n; jj++)
                        {
                            y = Ur[jj][j];
                            z = Ur[jj][i];
                            Ur[jj][j] = y * c + z * e;
                            Ur[jj][i] = z * c - y * e;
                        }
                    }

                    rv1[l] = 0.0f;
                    rv1[k] = f;
                    Sr[k] = x;
                }
            }

            // sort singular values descending and permute U, V columns accordingly
            for (i = 0; i < m - 1; i++)
            {
                int maxIdx = i;
                double maxVal = Sr[i];
                for (j = i + 1; j < m; j++)
                {
                    if (Sr[j] > maxVal)
                    {
                        maxVal = Sr[j];
                        maxIdx = j;
                    }
                }
                if (maxIdx != i)
                {
                    // swap S
                    var tS = Sr[i]; Sr[i] = Sr[maxIdx]; Sr[maxIdx] = tS;
                    // swap columns in U (n x m)
                    SwapColumns(Ur, i, maxIdx);
                    // swap columns in V (m x m)
                    SwapColumns(Vr, i, maxIdx);
                }
            }
            // Orthogonal factors do not depend on a positive common scale.
            this.Ur = Jagged.Zero(n, m);
            this.Vr = Jagged.Zero(m, m);
            this.Sr = new float[m];
            for (i = 0; i < n; i++)
                for (j = 0; j < m; j++) this.Ur[i][j] = (float)Ur[i][j];
            for (i = 0; i < m; i++)
            {
                this.Sr[i] = (float)(Sr[i] * inputScale);
                for (j = 0; j < m; j++) this.Vr[i][j] = (float)Vr[i][j];
            }
        }
        /// <summary>
        /// Copies a finite matrix and scales its largest magnitude to one before bidiagonalization.
        /// </summary>
        /// <param name="matrix">Real input matrix, with at least as many rows as columns.</param>
        /// <param name="scale">Receives the original maximum magnitude, or one for a zero matrix.</param>
        /// <returns>A scaled jagged copy; the caller must multiply the resulting singular values by scale.</returns>
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
        /// Swaps two columns in a jagged matrix <paramref name="M"/> (double[rows][cols]).
        /// No operation is performed if <paramref name="c1"/> equals <paramref name="c2"/>.
        /// </summary>
        /// <param name="M">Matrix represented as an array of row arrays (double[rows][cols])</param>
        /// <param name="c1">Index of the first column</param>
        /// <param name="c2">Index of the second column</param>
        private static void SwapColumns(double[][] M, int c1, int c2)
        {
            if (c1 == c2) return;
            int rows = M.Length;

            for (int r = 0; r < rows; r++)
            {
                var tmp = M[r][c1];
                M[r][c1] = M[r][c2];
                M[r][c2] = tmp;
            }
        }
        /// <summary>
        /// Computes sqrt(a*a + b*b) by a ratio to avoid squaring the larger magnitude.
        /// </summary>
        /// <param name="a">First finite real component.</param>
        /// <param name="b">Second finite real component.</param>
        /// <returns>The nonnegative Euclidean length, including zero for two zero components.</returns>
        private static double Hypotenuse(double a, double b)
        {
            a = Math.Abs(a); b = Math.Abs(b);
            if (a < b) { double temporary = a; a = b; b = temporary; }
            if (a == 0) return 0;
            double ratio = b / a;
            return a * Math.Sqrt(1 + ratio * ratio);
        }

        /// <summary>
        /// Transfers an algebraic sign to a magnitude for stable Householder and QR shifts.
        /// </summary>
        /// <param name="magnitude">Real magnitude donor.</param>
        /// <param name="sign">Sign donor; zero selects the nonnegative sign.</param>
        /// <returns>The absolute magnitude with the selected sign.</returns>
        private static double CopySign(double magnitude, double sign)
        {
            return sign < 0 ? -Math.Abs(magnitude) : Math.Abs(magnitude);
        }
        #endregion
    }
}
