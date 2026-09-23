using System;
using C = System.Numerics.Complex;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines eigenvalue decomposition.
    /// </summary>
    /// <remarks>
    /// Computes right eigenvectors satisfying A * V = V * D. For complex inputs, D is diagonal.
    /// For real inputs with complex-conjugate eigenvalues, V stores the real and imaginary parts
    /// in adjacent columns and D contains the corresponding real 2-by-2 blocks.
    /// Defective matrices need not have a complete set of independent eigenvectors.
    /// More information can be found on the website:
    /// <see href="https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix"/>.
    /// </remarks>
    public static class EVD
    {
        /// <summary>Computes the real EVD decomposition without modifying the inputs.</summary>
        /// <param name="matrix">Finite nonempty real square matrix.</param>
        /// <param name="eps">Relative convergence tolerance, clamped to [0,1] with a floor of eight double-precision rounding units.</param>
        /// <returns>Right eigenvectors V and eigenvalues D. Symmetric inputs have orthogonal V and ascending real D; other inputs use adjacent real columns for complex-conjugate eigenvectors.</returns>
        /// <exception cref="InvalidOperationException">The QR or QL iteration limit is reached before convergence.</exception>
        public static (float[,] V, Complex32[] D) Decompose(float[,] matrix, float eps = 1e-16f)
        {
            InternalMatrixMath.Validate(matrix, true);
            if (float.IsNaN(eps)) throw new ArgumentOutOfRangeException(nameof(eps));
            var d = Factor(matrix, eps);
            var values = new Complex32[d.Re.Length];
            for (int i = 0; i < values.Length; i++) values[i] = new Complex32((float)d.Re[i], (float)d.Im[i]);
            return (InternalMatrixMath.Real(d.V), values);
        }

        /// <summary>Computes right eigenvectors and eigenvalues of a complex square matrix.</summary>
        /// <param name="matrix">Finite nonempty square matrix, not modified.</param>
        /// <param name="eps">Relative convergence tolerance, clamped to [0,1] with a floor of eight double-precision rounding units.</param>
        /// <returns>V and D satisfying A V = V diag(D). Exactly Hermitian inputs have unitary V and ascending real D. Defective inputs need not have independent eigenvectors.</returns>
        public static (Complex32[,] V, Complex32[] D) Decompose(Complex32[,] matrix, float eps = 1e-16f)
        {
            if (float.IsNaN(eps)) throw new ArgumentOutOfRangeException(nameof(eps));
            var a = InternalMatrixMath.Copy(matrix, true);
            // Automatic dispatch must not erase a small imaginary eigenvalue by treating
            // a nearly Hermitian matrix as exactly Hermitian.
            if (InternalMatrixMath.IsHermitian(a, 0))
            {
                double tolerance = Math.Max(8 * InternalMatrixMath.Roundoff, Math.Min(1, Math.Max(0, eps)));
                var d = HermitianFactor(a, tolerance);
                var eigenvalues = new Complex32[d.D.Length];
                for (int i = 0; i < eigenvalues.Length; i++) eigenvalues[i] = new Complex32((float)d.D[i], 0);
                return (InternalMatrixMath.Single(d.V), eigenvalues);
            }
            var schur = Schur.Factor(a, eps);
            int n = a.GetLength(0);
            var alpha = new C[n];
            var beta = new double[n];
            var values = new Complex32[n];
            for (int i = 0; i < n; i++)
            {
                alpha[i] = schur.T[i, i]; beta[i] = 1;
                values[i] = new Complex32((float)alpha[i].Real, (float)alpha[i].Imaginary);
            }
            var vectors = TriangularVectors(schur.T, InternalMatrixMath.Eye(n), schur.Q, alpha, beta);
            return (InternalMatrixMath.Single(vectors), values);
        }

        /// <summary>Builds the complex diagonal eigenvalue matrix from existing eigenvalues.</summary>
        /// <param name="values">Eigenvalues in the order of the eigenvector columns.</param>
        /// <returns>diag(values), without repeating the eigenvalue calculation.</returns>
        public static Complex32[,] EigenvalueMatrix(Complex32[] values)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            return values.Diag();
        }

        /// <summary>Builds the real block eigenvalue matrix used with real eigenvector storage.</summary>
        /// <param name="values">Real eigenvalues or adjacent conjugate pairs, positive imaginary part first.</param>
        /// <returns>One-dimensional real blocks and two-dimensional blocks [a,b;-b,a].</returns>
        public static float[,] RealEigenvalueMatrix(Complex32[] values)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            int n = values.Length;
            var r = new float[n, n];
            for (int i = 0; i < n; i++)
            {
                r[i, i] = values[i].Real;
                if (values[i].Imag == 0) continue;
                if (values[i].Imag < 0 || i + 1 >= n || values[i + 1].Imag >= 0 ||
                    C.Abs(new C(values[i + 1].Real - (double)values[i].Real, values[i + 1].Imag + (double)values[i].Imag)) >
                        16 * InternalMatrixMath.SingleRoundoff * C.Abs(new C(values[i].Real, values[i].Imag)))
                    throw new ArgumentException("Complex eigenvalues must be adjacent conjugate pairs, positive imaginary part first.", nameof(values));
                r[i, i + 1] = values[i].Imag;
                // Homogeneous alpha/beta quotients can differ by a few rounding units within a conjugate pair.
                r[i + 1, i] = values[i + 1].Imag;
                r[i + 1, i + 1] = values[i + 1].Real;
                i++;
            }
            return r;
        }

        /// <summary>Computes the Hessenberg form separately when intermediate reduction data are required.</summary>
        /// <param name="matrix">Finite nonempty real square matrix.</param>
        /// <returns>The Hessenberg matrix of an independent reduction.</returns>
        public static float[,] HessenbergForm(float[,] matrix) => Hessenberg.Decompose(matrix).H;

        /// <summary>Computes the complex Hessenberg form separately from eigenvector calculation.</summary>
        /// <param name="matrix">Finite nonempty complex square matrix.</param>
        /// <returns>The Hessenberg matrix of an independent reduction.</returns>
        public static Complex32[,] HessenbergForm(Complex32[,] matrix) => Hessenberg.Decompose(matrix).H;

        /// <summary>Back-substitutes homogeneous eigenvectors of an upper triangular matrix pencil.</summary>
        /// <param name="s">First upper triangular factor.</param>
        /// <param name="t">Second upper triangular factor.</param>
        /// <param name="z">Right unitary transformation to original coordinates.</param>
        /// <param name="alpha">Eigenvalue numerators in diagonal order.</param>
        /// <param name="beta">Real nonnegative eigenvalue denominators.</param>
        /// <returns>Right eigenvectors normalized to maximum component magnitude one; repeated roots may yield dependent columns.</returns>
        internal static C[,] TriangularVectors(C[,] s, C[,] t, C[,] z, C[] alpha, double[] beta)
        {
            int n = s.GetLength(0);
            var vectors = new C[n, n];
            for (int k = 0; k < n; k++)
            {
                // Normalize coefficients before multiplication; homogeneous pairs may span very different units.
                double pairScale = Math.Max(C.Abs(alpha[k]), Math.Abs(beta[k]));
                C numerator = pairScale == 0 ? C.Zero : alpha[k] / pairScale;
                double denominator = pairScale == 0 ? 0 : beta[k] / pairScale;
                var y = new C[n]; y[k] = 1;
                for (int i = k - 1; i >= 0; i--)
                {
                    C sum = 0;
                    double local = C.Abs(denominator * s[i, i]) + C.Abs(numerator * t[i, i]);
                    for (int j = i + 1; j <= k; j++)
                    {
                        C coefficient = denominator * s[i, j] - numerator * t[i, j];
                        sum += coefficient * y[j];
                        if (y[j] != C.Zero) local = Math.Max(local, C.Abs(coefficient));
                    }
                    // Only the active row determines a near-repeated-root safeguard.
                    // A large eigenvalue in an independent block must not perturb this equation.
                    double floor = Math.Max(1e-300, 16 * InternalMatrixMath.Roundoff * local);
                    C diagonal = denominator * s[i, i] - numerator * t[i, i];
                    if (C.Abs(diagonal) < floor)
                    {
                        if (sum == C.Zero) { y[i] = 0; continue; }
                        diagonal = C.Abs(diagonal) == 0 ? new C(floor, 0) : diagonal / C.Abs(diagonal) * floor;
                    }
                    y[i] = -sum / diagonal;
                    double magnitude = C.Abs(y[i]);
                    if (magnitude > 1e100)
                        for (int j = i; j <= k; j++) y[j] /= magnitude;
                }
                double largest = 0;
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j <= k; j++) vectors[i, k] += z[i, j] * y[j];
                    largest = Math.Max(largest, C.Abs(vectors[i, k]));
                }
                if (largest > 0)
                    for (int i = 0; i < n; i++) vectors[i, k] /= largest;
            }
            return vectors;
        }

        /// <summary>Factors a real input with symmetric QL or nonsymmetric Hessenberg/QR iteration.</summary>
        /// <remarks>All work buffers belong to this call; conjugate eigenvectors use adjacent real columns.</remarks>
        private static (double[][] V, double[] Re, double[] Im) Factor(float[,] a, double eps)
        {
            int n = a.GetLength(0);
            var real = new double[n];
            var imaginary = new double[n];
            eps = Math.Max(8 * InternalMatrixMath.Roundoff, Math.Min(1, Math.Max(0, eps)));
            double[][] vectors;
            if (Matrice.IsSymmetric(a))
                vectors = SymmetricFactor(a, real, imaginary, eps);
            else
            {
                vectors = InternalMatrixMath.CreateJagged(n, n);
                var hessenberg = InternalMatrixMath.CopyJagged(a);
                ReduceToHessenberg(hessenberg, vectors);
                ReduceHessenbergToEigenvectors(hessenberg, vectors, real, imaginary, eps);
            }
            return (vectors, real, imaginary);
        }

        /// <summary>Uses Householder tridiagonalization and shared QL iteration for symmetric real inputs.</summary>
        private static double[][] SymmetricFactor(float[,] a, double[] real, double[] imaginary, double eps)
        {
            int n = a.GetLength(0);
            var components = n > 1
                ? InternalMatrixMath.ConnectedComponents(n, (i, j) => a[i, j] != 0)
                : Array.Empty<int[]>();
            if (components.Length > 1)
            {
                var vectors = InternalMatrixMath.CreateJagged(n, n);
                int start = 0;
                foreach (var indices in components)
                {
                    int size = indices.Length;
                    var input = new float[size, size];
                    for (int i = 0; i < size; i++)
                        for (int j = 0; j < size; j++) input[i, j] = a[indices[i], indices[j]];
                    var part = Factor(input, eps);
                    for (int j = 0; j < size; j++)
                    {
                        real[start + j] = part.Re[j];
                        for (int i = 0; i < size; i++) vectors[indices[i]][start + j] = part.V[i][j];
                    }
                    start += size;
                }
                SortSymmetricEigenpairs(real, (i, j) => InternalMatrixMath.SwapColumns(vectors, i, j));
                return vectors;
            }
            var matrices = InternalMatrixMath.CopyJagged(a);
            Tridiagonalize(matrices, real, imaginary);
            DiagonalizeTridiagonal(real, imaginary, eps,
                (i, j, c, t) => InternalMatrixMath.RotateColumns(matrices, i, j, c, t),
                (i, j) => InternalMatrixMath.SwapColumns(matrices, i, j));
            return matrices;
        }

        /// <summary>Uses the real tridiagonal/QL strategy for exactly Hermitian complex inputs.</summary>
        private static (C[,] V, double[] D) HermitianFactor(C[,] a, double eps)
        {
            int n = a.GetLength(0);
            var components = InternalMatrixMath.ConnectedComponents(n, (i, j) => a[i, j] != C.Zero);
            if (components.Length > 1)
            {
                var vectors = new C[n, n];
                var values = new double[n];
                int column = 0;
                foreach (var indices in components)
                {
                    var block = new C[indices.Length, indices.Length];
                    for (int i = 0; i < indices.Length; i++)
                        for (int j = 0; j < indices.Length; j++) block[i, j] = a[indices[i], indices[j]];
                    var part = HermitianFactor(block, eps);
                    for (int j = 0; j < indices.Length; j++)
                    {
                        values[column + j] = part.D[j];
                        for (int i = 0; i < indices.Length; i++) vectors[indices[i], column + j] = part.V[i, j];
                    }
                    column += indices.Length;
                }
                SortSymmetricEigenpairs(values, (i, j) => InternalMatrixMath.SwapColumns(vectors, i, j));
                return (vectors, values);
            }
            var reduction = Householder.Tridiagonalize(a);
            var d = new double[n];
            var e = new double[n];
            C phase = C.One;
            for (int j = 0; j < n; j++)
            {
                d[j] = reduction.H[j, j].Real;
                if (j > 0)
                {
                    e[j] = C.Abs(reduction.H[j, j - 1]);
                    phase *= InternalMatrixMath.Phase(reduction.H[j, j - 1]);
                }
                for (int i = 0; i < n; i++) reduction.P[i, j] *= phase;
            }
            DiagonalizeTridiagonal(d, e, eps,
                (i, j, c, sine) => InternalMatrixMath.RotateColumns(reduction.P, i, j, c, sine),
                (i, j) => InternalMatrixMath.SwapColumns(reduction.P, i, j));
            return (reduction.P, d);
        }

        /// <summary>Shared real tridiagonal QL iteration derived from EISPACK tql2.</summary>
        /// <remarks>
        /// Derived from tql2 by Bowdler, Martin, Reinsch, and Wilkinson, Handbook for Auto. Comp., Vol. II, Linear Algebra.
        /// Im[0] is unused on entry; Im[i] is the subdiagonal at row i. Transformations are accumulated in the original scalar domain.
        /// </remarks>
        private static void DiagonalizeTridiagonal(double[] Re, double[] Im, double eps,
            Action<int, int, double, double> rotate, Action<int, int> swap)
        {
            int n = Re.Length;

            double f = 0;
            double tst1 = 0;
            int blockEnd = -1;
            int i, l, iter, m;
            double g, p, r, dl1, h;
            double c, c2, c3, el1, s, s2;

            for (i = 1; i < n; i++)
                Im[i - 1] = Im[i];

            Im[n - 1] = 0;

            for (l = 0; l < n; l++)
            {
                // Independent tridiagonal blocks need independent shifts and scales.
                // Carrying either across an exact zero coupling can erase a small spectrum.
                if (l > blockEnd)
                {
                    f = tst1 = 0;
                    blockEnd = l;
                    while (blockEnd < n - 1 && Im[blockEnd] != 0) blockEnd++;
                }
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
                        if (++iter > 1000) throw new InvalidOperationException("Symmetric/Hermitian EVD failed to converge.");

                        // Compute implicit shift
                        g = Re[l];
                        p = (Re[l + 1] - g) / (2 * Im[l]);
                        r = InternalMatrixMath.Hypotenuse(p, 1);
                        if (p < 0)
                        {
                            r = -r;
                        }

                        Re[l] = Im[l] / (p + r);
                        Re[l + 1] = Im[l] * (p + r);
                        dl1 = Re[l + 1];
                        h = g - Re[l];
                        for (i = l + 2; i <= blockEnd; i++)
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
                            r = InternalMatrixMath.Hypotenuse(p, Im[i]);
                            Im[i + 1] = s * r;
                            s = Im[i] / r;
                            c = p / r;
                            p = c * Re[i] - s * g;
                            Re[i + 1] = h + s * (c * g + s * Re[i]);

                            // Accumulate transformation.
                            rotate(i, i + 1, c, -s);
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

            SortSymmetricEigenpairs(Re, swap);
        }

        /// <summary>Orders a symmetric or Hermitian spectrum and its eigenvectors together.</summary>
        private static void SortSymmetricEigenpairs(double[] values, Action<int, int> swap)
        {
            for (int i = 0; i < values.Length - 1; i++)
            {
                int best = i;
                for (int j = i + 1; j < values.Length; j++) if (values[j] < values[best]) best = j;
                if (best == i) continue;
                double value = values[i]; values[i] = values[best]; values[best] = value;
                swap(i, best);
            }
        }

        /// <summary>
        /// Symmetric Householder reduction to tridiagonal form.
        /// This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch, and Wilkinson,
        /// Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
        /// </summary>
        private static void Tridiagonalize(double[][] matrices, double[] Re, double[] Im)
        {
            int n = matrices.Length;
            int i, j, k;

            for (j = 0; j < n; j++)
            {
                Re[j] = matrices[n - 1][j];
            }

            double scale, h, f, g, hh;

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
                    g = System.Math.Sqrt(h);
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
        /// Nonsymmetric reduction to Hessenberg form.
        /// This is derived from the Algol procedures orthes and ortran, by Martin and Wilkinson,
        /// Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutines in EISPACK.
        /// </summary>
        private static void ReduceToHessenberg(double[][] hessenberg, double[][] matrices)
        {
            int n = hessenberg.Length;
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

                        // Divide in two stages to avoid underflow in the denominator product.
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
        private static void ReduceHessenbergToEigenvectors(double[][] hessenberg, double[][] matrices, double[] Re, double[] Im, double eps)
        {
            int nn = hessenberg.Length;
            int n = nn - 1;
            int low = 0;
            int high = nn - 1;
            double exshift = 0;
            double p = 0;
            double q = 0;
            double r = 0;
            double s = 0;
            double z = 0;
            double t;
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

                    if (double.IsNaN(s))
                        break;

                    if (System.Math.Abs(hessenberg[l][l - 1]) <= eps * s)
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
                        x = y = 0.75 * s;
                        w = (-0.4375) * s * s;
                    }

                    // MATLAB's new ad hoc shift
                    if (iter == 30)
                    {
                        s = (y - x) / 2;
                        s = s * s + w;
                        if (s > 0)
                        {
                            s = System.Math.Sqrt(s);
                            if (y < x) s = -s;
                            s = x - w / ((y - x) / 2 + s);
                            for (i = low; i <= n; i++)
                                hessenberg[i][i] -= s;
                            exshift += s;
                            x = y = w = 0.964;
                        }
                    }

                    if (++iter > 1000) throw new InvalidOperationException("Real nonsymmetric EVD failed to converge.");

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

                    // Double QR step involving rows l:n and columns m:n.
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

                        s = System.Math.Sqrt(p * p + q * q + r * r);
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
                        InternalMatrixMath.DivideComplex(0, -hessenberg[n - 1][n], hessenberg[n - 1][n - 1] - p, q, ref hessenberg[n - 1][n - 1], ref hessenberg[n - 1][n]);
                    }

                    hessenberg[n][n - 1] = 0;
                    hessenberg[n][n] = 1;
                    for (i = n - 2; i >= 0; i--)
                    {
                        double ra, sa, vr, vi;
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
                                InternalMatrixMath.DivideComplex(-ra, -sa, w, q, ref hessenberg[i][n - 1], ref hessenberg[i][n]);
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
                                InternalMatrixMath.DivideComplex(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi, ref hessenberg[i][n - 1], ref hessenberg[i][n]);
                                if (System.Math.Abs(x) > (System.Math.Abs(z) + System.Math.Abs(q)))
                                {
                                    hessenberg[i + 1][n - 1] = (-ra - w * hessenberg[i][n - 1] + q * hessenberg[i][n]) / x;
                                    hessenberg[i + 1][n] = (-sa - w * hessenberg[i][n] - q * hessenberg[i][n - 1]) / x;
                                }
                                else
                                {
                                    InternalMatrixMath.DivideComplex(-r - y * hessenberg[i][n - 1], -s - y * hessenberg[i][n], z, q, ref hessenberg[i + 1][n - 1], ref hessenberg[i + 1][n]);
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
    }
}
