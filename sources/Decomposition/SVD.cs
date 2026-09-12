using System;
using C = System.Numerics.Complex;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines singular value decomposition.
    /// </summary>
    /// <remarks>
    /// Represents a rectangular matrix as A = U * diag(S) * V^H for complex inputs,
    /// or A = U * diag(S) * V^T for real inputs. S contains nonnegative singular values.
    /// More information can be found on the website:
    /// <see href="https://en.wikipedia.org/wiki/Singular_value_decomposition"/>.
    /// </remarks>
    public static class SVD
    {
        /// <summary>Computes the economy real SVD, A = U diag(S) V^T, without modifying the input.</summary>
        /// <param name="matrix">Finite nonempty m by n matrix.</param>
        /// <param name="iterations">Positive maximum QR sweeps per singular value.</param>
        /// <returns>U of size m by k, descending nonnegative S of length k, and V of size n by k, where k=min(m,n).</returns>
        public static (float[,] U, float[] S, float[,] V) Decompose(float[,] matrix, int iterations = 10)
        {
            var d = Factor(InternalMatrixMath.CopyReal(matrix), iterations);
            return (InternalMatrixMath.Real(d.U), InternalMatrixMath.Single(d.S), InternalMatrixMath.Real(d.V));
        }

        /// <summary>Computes the economy complex SVD, A = U diag(S) V^H, using one-sided Jacobi sweeps.</summary>
        /// <param name="matrix">Finite nonempty m by n complex matrix, not modified.</param>
        /// <param name="iterations">Positive maximum number of cyclic Jacobi sweeps.</param>
        /// <returns>U of size m by k, descending nonnegative S of length k, and V of size n by k, where k=min(m,n).</returns>
        /// <exception cref="InvalidOperationException">The iteration limit is reached before convergence.</exception>
        public static (Complex32[,] U, float[] S, Complex32[,] V) Decompose(Complex32[,] matrix, int iterations = 50)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix), iterations);
            return (InternalMatrixMath.Single(d.U), InternalMatrixMath.Single(d.S), InternalMatrixMath.Single(d.V));
        }

        /// <summary>Constructs the Moore-Penrose inverse from existing real economy SVD factors.</summary>
        /// <param name="u">Left singular vectors, m by k.</param>
        /// <param name="s">Nonnegative singular values of length k.</param>
        /// <param name="v">Right singular vectors, n by k.</param>
        /// <param name="tolerance">Relative rank cutoff; -1 uses max(m,n)*2^-23.</param>
        /// <returns>The n by m pseudoinverse. No decomposition is repeated.</returns>
        public static float[,] PseudoInverse(float[,] u, float[] s, float[,] v, float tolerance = -1)
            => InternalMatrixMath.Real(Inverse(InternalMatrixMath.Copy(u), s, InternalMatrixMath.Copy(v), tolerance));

        /// <summary>Constructs the Moore-Penrose inverse from existing complex economy SVD factors.</summary>
        /// <param name="u">Left singular vectors, m by k.</param>
        /// <param name="s">Nonnegative singular values of length k.</param>
        /// <param name="v">Right singular vectors, n by k.</param>
        /// <param name="tolerance">Relative rank cutoff; -1 uses max(m,n)*2^-23.</param>
        /// <returns>V diag(S+) U^H, of size n by m. No decomposition is repeated.</returns>
        public static Complex32[,] PseudoInverse(Complex32[,] u, float[] s, Complex32[,] v, float tolerance = -1)
            => InternalMatrixMath.Single(Inverse(InternalMatrixMath.Copy(u), s, InternalMatrixMath.Copy(v), tolerance));

        /// <summary>Counts singular values exceeding a relative cutoff.</summary>
        /// <param name="s">Finite nonnegative singular values.</param>
        /// <param name="tolerance">Relative cutoff; -1 uses length(S)*2^-23. Supply max(m,n)*2^-23 for rectangular inputs.</param>
        /// <returns>The numerical rank.</returns>
        public static int Rank(float[] s, float tolerance = -1)
        {
            double largest = ValidateValues(s);
            double cutoff = Cutoff(tolerance, s.Length) * largest;
            int rank = 0;
            foreach (float value in s) if (value > cutoff) rank++;
            return rank;
        }

        /// <summary>Returns the spectral norm from existing singular values.</summary>
        /// <param name="s">Finite nonnegative singular values.</param>
        /// <returns>The largest singular value, or zero for an empty sequence.</returns>
        public static float Norm(float[] s) => (float)ValidateValues(s);

        /// <summary>Returns the spectral condition number from existing singular values.</summary>
        /// <param name="s">Finite nonnegative singular values.</param>
        /// <returns>max(S)/min(S), or positive infinity if the sequence is empty or contains zero.</returns>
        public static float ConditionNumber(float[] s)
        {
            double largest = ValidateValues(s), smallest = double.PositiveInfinity;
            foreach (float value in s) smallest = Math.Min(smallest, value);
            return s.Length == 0 || smallest == 0 ? float.PositiveInfinity : (float)(largest / smallest);
        }

        /// <summary>Validates a relative rank cutoff and selects the single-precision default.</summary>
        /// <param name="tolerance">Nonnegative relative cutoff, or -1 for the default.</param>
        /// <param name="dimension">Dimension used for the default error scale.</param>
        /// <returns>A nonnegative relative cutoff.</returns>
        private static double Cutoff(float tolerance, int dimension)
        {
            if (float.IsNaN(tolerance) || float.IsInfinity(tolerance) || (tolerance < 0 && tolerance != -1))
                throw new ArgumentOutOfRangeException(nameof(tolerance));
            return tolerance == -1 ? dimension * InternalMatrixMath.SingleRoundoff : tolerance;
        }

        /// <summary>Validates finite nonnegative singular values without requiring sorted order.</summary>
        /// <param name="s">Singular-value vector.</param>
        /// <returns>The largest value, or zero for an empty vector.</returns>
        private static double ValidateValues(float[] s)
        {
            if (s == null) throw new ArgumentNullException(nameof(s));
            double largest = 0;
            foreach (float value in s)
            {
                if (value < 0 || float.IsNaN(value) || float.IsInfinity(value))
                    throw new ArgumentException("Singular values must be finite and nonnegative.", nameof(s));
                largest = Math.Max(largest, value);
            }
            return largest;
        }

        /// <summary>Accumulates a truncated pseudoinverse with conjugation and double-precision division.</summary>
        /// <param name="u">Left factors.</param>
        /// <param name="s">Singular values.</param>
        /// <param name="v">Right factors.</param>
        /// <param name="tolerance">Relative cutoff or -1 for the dimension-dependent default.</param>
        /// <returns>The pseudoinverse; inputs are not modified.</returns>
        private static C[,] Inverse(C[,] u, float[] s, C[,] v, float tolerance)
        {
            double largest = ValidateValues(s);
            if (u.GetLength(1) != s.Length || v.GetLength(1) != s.Length)
                throw new ArgumentException("The factor dimensions must agree with the number of singular values.");
            int m = u.GetLength(0), n = v.GetLength(0);
            double threshold = Cutoff(tolerance, Math.Max(m, n)) * largest;
            var result = new C[n, m];
            for (int k = 0; k < s.Length; k++)
                if (s[k] > threshold)
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < m; j++)
                            result[i, j] += (v[i, k] / s[k]) * C.Conjugate(u[j, k]);
            return result;
        }

        /// <summary>Orthogonalizes complex columns using unitary plane rotations without forming A^H A.</summary>
        /// <param name="a">Private work matrix.</param>
        /// <param name="iterations">Maximum positive sweep count.</param>
        /// <returns>Double-precision economy singular factors and descending singular values.</returns>
        internal static (C[,] U, double[] S, C[,] V) Factor(C[,] a, int iterations = 50)
        {
            if (iterations < 1) throw new ArgumentOutOfRangeException(nameof(iterations));
            int m = a.GetLength(0), n = a.GetLength(1);
            if (m < n)
            {
                var wide = Factor(InternalMatrixMath.Adjoint(a), iterations);
                return (wide.V, wide.S, wide.U);
            }
            var v = InternalMatrixMath.Eye(n);
            bool converged = n < 2;
            for (int sweep = 0; sweep < iterations && !converged; sweep++)
            {
                converged = true;
                for (int p = 0; p < n - 1; p++)
                    for (int q = p + 1; q < n; q++)
                    {
                        double np = InternalMatrixMath.ColumnNorm(a, p), nq = InternalMatrixMath.ColumnNorm(a, q);
                        double scale = Math.Max(np, nq);
                        if (np == 0 || nq == 0) continue;
                        double ap = np / scale, aq = nq / scale;
                        C dot = 0;
                        for (int i = 0; i < m; i++) dot += C.Conjugate(a[i, p] / np) * (a[i, q] / nq);
                        double correlation = C.Abs(dot);
                        if (correlation <= 8 * InternalMatrixMath.Roundoff * Math.Max(1, m)) continue;
                        // Normalize the two-column Gram entries locally, preserving isolated tiny singular values.
                        double cross = correlation * ap * aq;
                        if (cross == 0) continue;
                        double tau = (aq * aq - ap * ap) / (2 * cross);
                        double t = (tau >= 0 ? 1 : -1) / (Math.Abs(tau) + Math.Sqrt(1 + tau * tau));
                        if (t == 0) continue;
                        double c = 1 / Math.Sqrt(1 + t * t), sn = t * c;
                        C phase = dot / correlation;
                        // Express the Jacobi update in the shared complex rotation convention.
                        C sine = -sn * C.Conjugate(phase);
                        InternalMatrixMath.RotateColumns(a, p, q, c, sine);
                        InternalMatrixMath.RotateColumns(v, p, q, c, sine);
                        converged = false;
                    }
            }
            if (!converged) throw new InvalidOperationException("Complex SVD failed to converge within the Jacobi sweep limit.");
            var singular = new double[n];
            for (int j = 0; j < n; j++) singular[j] = InternalMatrixMath.ColumnNorm(a, j);
            for (int p = 0; p < n; p++)
            {
                int best = p;
                for (int q = p + 1; q < n; q++) if (singular[q] > singular[best]) best = q;
                if (best == p) continue;
                double d = singular[p]; singular[p] = singular[best]; singular[best] = d;
                InternalMatrixMath.SwapColumns(a, p, best);
                InternalMatrixMath.SwapColumns(v, p, best);
            }
            var u = new C[m, n];
            for (int j = 0; j < n; j++)
            {
                if (singular[j] > 0)
                    for (int i = 0; i < m; i++) u[i, j] = a[i, j] / singular[j];
                else
                {
                    var column = InternalMatrixMath.Complete(u, j);
                    InternalMatrixMath.SetColumn(u, j, column);
                }
            }
            return (u, singular, v);
        }

        /// <summary>Consumes a private real buffer without narrowing intermediate singular factors.</summary>
        internal static (double[][] U, double[] S, double[][] V) Factor(double[][] a, int iterations = 50)
        {
            var work = new RealWorkspace(a, iterations);
            return (work.U, work.S, work.V);
        }

        /// <summary>Owns the real algorithm work buffers for one call only.</summary>
        private sealed class RealWorkspace
        {
            #region Private data
            private int n, m;
            private int iterations;
            private double[][] Ur;
            private double[][] Vr;
            private double[] Sr;
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
            public RealWorkspace(double[][] A, int iterations = 10)
            {
                if (A == null) throw new ArgumentNullException(nameof(A));
                if (A.Length == 0 || A[0].Length == 0)
                    throw new ArgumentException("The matrix must be nonempty.", nameof(A));
                if (iterations < 1) throw new ArgumentOutOfRangeException(nameof(iterations), "The iteration limit must be positive.");
                // set:
                this.iterations = iterations;
                this.n = A.Length;
                this.m = A[0].Length;

                // options:
                if (n < m)
                {
                    this.reversed = true;
                    this.n = A[0].Length;
                    this.m = A.Length;
                    this.svdcmp(InternalMatrixMath.Transpose(A));
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
            public double[][] U
            {
                get
                {
                    return reversed ? Vr : Ur;
                }
            }
            /// <summary>
            /// Gets singular values.
            /// </summary>
            public double[] S
            {
                get { return Sr; }
            }
            /// <summary>
            /// Gets the right vectors.
            /// </summary>
            public double[][] V
            {
                get
                {
                    return reversed ? Ur : Vr;
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
            /// beforehand if needed). The method works on an internal copy (jagged buffers).
            /// </param>
            /// <remarks>
            /// Uses jagged arrays for speed. Columns of U and V are orthonormal. The number
            /// of QR sweeps is limited by the instance field <see cref="iterations"/>.
            /// </remarks>
            private void svdcmp(double[][] A)
            {
                double inputScale = InternalMatrixMath.Max(A);
                if (inputScale == 0) inputScale = 1;
                InternalMatrixMath.Divide(A, inputScale);
                var Ur = A;
                var Sr = new double[m];
                var Vr = InternalMatrixMath.CreateJagged(m, m);
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
                            g = -InternalMatrixMath.CopySign(Math.Sqrt(e), f);
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
                            g = -InternalMatrixMath.CopySign(Math.Sqrt(e), f);
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
                                h = InternalMatrixMath.Hypotenuse(f, g);
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
                        g = InternalMatrixMath.Hypotenuse(f, 1.0f);
                        f = ((x - z) * (x + z) + h * ((y / (f + InternalMatrixMath.CopySign(g, f))) - h)) / x;

                        // next QR transformation
                        c = e = 1.0f;

                        for (j = l; j <= nm; j++)
                        {
                            i = j + 1;
                            g = rv1[i];
                            y = Sr[i];
                            h = e * g;
                            g = c * g;
                            z = InternalMatrixMath.Hypotenuse(f, h);
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

                            z = InternalMatrixMath.Hypotenuse(f, h);
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
                        InternalMatrixMath.SwapColumns(Ur, i, maxIdx);
                        // swap columns in V (m x m)
                        InternalMatrixMath.SwapColumns(Vr, i, maxIdx);
                    }
                }
                // Orthogonal factors do not depend on a positive common scale.
                this.Ur = Ur;
                this.Vr = Vr;
                this.Sr = Sr;
                for (i = 0; i < m; i++) Sr[i] *= inputScale;
            }

            #endregion

        }
    }
}
