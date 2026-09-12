using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides real and complex generalized Schur decomposition.</summary>
    public static class QZ
    {
        /// <summary>Computes the real generalized Schur factors A = Q S Z^T and B = Q T Z^T.</summary>
        /// <param name="a">Finite nonempty square matrix A.</param>
        /// <param name="b">Finite square matrix B of the same order.</param>
        /// <param name="eps">Relative deflation tolerance with a roundoff floor.</param>
        /// <returns>Orthogonal Q and Z, quasi-triangular S, and upper triangular T.</returns>
        public static (float[,] Q, float[,] S, float[,] T, float[,] Z) Decompose(float[,] a, float[,] b, float eps = 1e-16f)
        {
            var s = InternalRealMatrixMath.Copy(a, true);
            var t = InternalRealMatrixMath.Copy(b, true);
            if (a.GetLength(0) != b.GetLength(0)) throw new ArgumentException("The matrices must have equal orders.");
            if (float.IsNaN(eps)) throw new ArgumentOutOfRangeException(nameof(eps));
            int n = a.GetLength(0), error = 0;
            // Store Q^T so each left reflection updates contiguous rows.
            var q = InternalRealMatrixMath.Eye(n);
            var z = InternalRealMatrixMath.Eye(n);
            GEVD.ReduceRealPencil(s, t, eps, q, z, ref error);
            if (error != 0) throw new InvalidOperationException("Real QZ decomposition failed to converge.");
            return (InternalRealMatrixMath.RealTranspose(q), InternalRealMatrixMath.Real(s),
                    InternalRealMatrixMath.Real(t), InternalRealMatrixMath.Real(z));
        }

        /// <summary>Computes the complex generalized Schur factors A = Q S Z^H and B = Q T Z^H.</summary>
        /// <param name="a">Finite nonempty square matrix A.</param>
        /// <param name="b">Finite square matrix B of the same order; it may be singular.</param>
        /// <param name="eps">Relative deflation tolerance with a roundoff floor.</param>
        /// <param name="iterations">Positive maximum QZ steps between successive deflations.</param>
        /// <returns>Unitary Q and Z and upper triangular S and T; T has real nonnegative diagonal entries.</returns>
        public static (Complex32[,] Q, Complex32[,] S, Complex32[,] T, Complex32[,] Z) Decompose(
            Complex32[,] a, Complex32[,] b, float eps = 1e-16f, int iterations = 1000)
        {
            if (float.IsNaN(eps)) throw new ArgumentOutOfRangeException(nameof(eps));
            var d = Factor(InternalMatrixMath.Copy(a, true), InternalMatrixMath.Copy(b, true), eps, iterations);
            return (InternalMatrixMath.Single(d.Q), InternalMatrixMath.Single(d.S), InternalMatrixMath.Single(d.T), InternalMatrixMath.Single(d.Z));
        }

        /// <summary>Reduces a complex pencil by unitary Hessenberg-triangular reduction and implicit single-shift QZ.</summary>
        /// <param name="a">Private first square matrix.</param>
        /// <param name="b">Private second square matrix of the same order.</param>
        /// <param name="eps">Relative deflation tolerance.</param>
        /// <param name="iterations">Maximum steps between deflations.</param>
        /// <returns>Double-precision generalized Schur factors, including infinite and indeterminate diagonal pairs.</returns>
        /// <remarks>Zero diagonals of B are chased without inverting B. See the LAPACK ZHGEQZ algorithm and Moler-Stewart QZ method.</remarks>
        internal static (C[,] Q, C[,] S, C[,] T, C[,] Z) Factor(C[,] a, C[,] b, double eps = 1e-16, int iterations = 1000)
        {
            int n = a.GetLength(0);
            if (b.GetLength(0) != n) throw new ArgumentException("The matrices must have equal orders.");
            if (iterations < 1) throw new ArgumentOutOfRangeException(nameof(iterations));
            double scaleA = InternalMatrixMath.Max(a), scaleB = InternalMatrixMath.Max(b);
            if (scaleA == 0) scaleA = 1;
            if (scaleB == 0) scaleB = 1;
            InternalMatrixMath.Divide(a, scaleA);
            InternalMatrixMath.Divide(b, scaleB);
            var qr = QR.Factor(b);
            b = qr.R;
            var q = qr.Q;
            var z = InternalMatrixMath.Eye(n);
            a = InternalMatrixMath.Multiply(InternalMatrixMath.Adjoint(q), a);
            for (int k = 0; k < n - 2; k++)
                for (int i = n - 1; i > k + 1; i--)
                {
                    var left = InternalMatrixMath.Givens(a[i - 1, k], a[i, k]);
                    LeftPair(a, b, q, i - 1, i, left.C, left.S);
                    a[i, k] = 0;
                    var right = InternalMatrixMath.Givens(b[i, i], b[i, i - 1]);
                    RightPair(a, b, z, i, i - 1, right.C, right.S);
                    b[i, i - 1] = 0;
                }
            double tolerance = Math.Max(8 * InternalMatrixMath.Roundoff, Math.Min(1, Math.Max(0, eps)));
            double bTolerance = tolerance * InternalMatrixMath.Max(b);
            int high = n - 1, steps = 0;
            while (high >= 0)
            {
                if (high == 0 || SmallSubdiagonal(a, high, tolerance))
                {
                    if (high > 0) a[high, high - 1] = 0;
                    high--; steps = 0; continue;
                }
                if (++steps > iterations) throw new InvalidOperationException("Complex QZ decomposition failed to converge.");
                if (C.Abs(b[high, high]) <= bTolerance)
                {
                    b[high, high] = 0;
                    var r = InternalMatrixMath.Givens(a[high, high], a[high, high - 1]);
                    RightPair(a, b, z, high, high - 1, r.C, r.S);
                    a[high, high - 1] = 0;
                    high--; steps = 0; continue;
                }
                int low = 0;
                bool chased = false;
                for (int j = high - 1; j >= 0; j--)
                {
                    bool split = j == 0 || SmallSubdiagonal(a, j, tolerance);
                    if (split && j > 0) a[j, j - 1] = 0;
                    if (C.Abs(b[j, j]) <= bTolerance)
                    {
                        b[j, j] = 0;
                        if (split)
                        {
                            var r = InternalMatrixMath.Givens(a[j, j], a[j + 1, j]);
                            LeftPair(a, b, q, j, j + 1, r.C, r.S);
                            a[j + 1, j] = 0;
                        }
                        else
                        {
                            // Move a zero B diagonal to the trailing corner while removing each Hessenberg bulge.
                            for (int k = j; k < high; k++)
                            {
                                var left = InternalMatrixMath.Givens(b[k, k + 1], b[k + 1, k + 1]);
                                LeftPair(a, b, q, k, k + 1, left.C, left.S);
                                b[k + 1, k + 1] = 0;
                                var right = InternalMatrixMath.Givens(a[k + 1, k], a[k + 1, k - 1]);
                                RightPair(a, b, z, k, k - 1, right.C, right.S);
                                a[k + 1, k - 1] = 0;
                            }
                        }
                        chased = true; break;
                    }
                    if (split) { low = j; break; }
                }
                if (chased) continue;
                C u12 = b[high - 1, high] / b[high, high];
                C d11 = a[high - 1, high - 1] / b[high - 1, high - 1];
                C d21 = a[high, high - 1] / b[high - 1, high - 1];
                C d12 = a[high - 1, high] / b[high, high];
                C d22 = a[high, high] / b[high, high];
                C bottom = d22 - u12 * d21;
                C center = (d11 + bottom) / 2;
                C root = C.Sqrt(center * center + d12 * d21 - d11 * d22);
                C shift1 = center + root, shift2 = center - root;
                C shift = C.Abs(shift1 - bottom) < C.Abs(shift2 - bottom) ? shift1 : shift2;
                if (steps % 10 == 0) shift = bottom + new C(0.75, 0.25) * C.Abs(d21);
                var rotation = InternalMatrixMath.Givens(a[low, low] - shift * b[low, low], a[low + 1, low]);
                for (int j = low; j < high; j++)
                {
                    if (j > low) rotation = InternalMatrixMath.Givens(a[j, j - 1], a[j + 1, j - 1]);
                    LeftPair(a, b, q, j, j + 1, rotation.C, rotation.S);
                    if (j > low) a[j + 1, j - 1] = 0;
                    var right = InternalMatrixMath.Givens(b[j + 1, j + 1], b[j + 1, j]);
                    RightPair(a, b, z, j + 1, j, right.C, right.S);
                    b[j + 1, j] = 0;
                }
            }
            for (int j = 0; j < n; j++)
            {
                double magnitude = C.Abs(b[j, j]);
                C phase = magnitude == 0 ? C.One : C.Conjugate(b[j, j]) / magnitude;
                for (int i = 0; i < n; i++)
                {
                    a[i, j] = i > j ? C.Zero : a[i, j] * phase * scaleA;
                    b[i, j] = i > j ? C.Zero : b[i, j] * phase * scaleB;
                    z[i, j] *= phase;
                }
                b[j, j] = magnitude * scaleB;
            }
            return (q, a, b, z);
        }

        /// <summary>Tests a Hessenberg subdiagonal against a local relative scale.</summary>
        /// <param name="a">Hessenberg matrix.</param>
        /// <param name="i">Subdiagonal row, greater than zero.</param>
        /// <param name="tolerance">Relative deflation threshold.</param>
        /// <returns>Whether the entry is negligible on its local scale.</returns>
        private static bool SmallSubdiagonal(C[,] a, int i, double tolerance)
        {
            double scale = C.Abs(a[i - 1, i - 1]) + C.Abs(a[i, i]);
            if (scale == 0) scale = C.Abs(a[i - 1, i]) + (i > 1 ? C.Abs(a[i - 1, i - 2]) : 0);
            return C.Abs(a[i, i - 1]) <= tolerance * scale;
        }

        /// <summary>Applies a left plane rotation to both matrices and updates Q.</summary>
        /// <param name="a">First work matrix.</param>
        /// <param name="b">Second work matrix.</param>
        /// <param name="q">Left accumulator.</param>
        /// <param name="i">First row.</param>
        /// <param name="j">Second row.</param>
        /// <param name="c">Real cosine.</param>
        /// <param name="s">Complex sine.</param>
        private static void LeftPair(C[,] a, C[,] b, C[,] q, int i, int j, double c, C s)
        {
            InternalMatrixMath.RotateRows(a, i, j, c, s);
            InternalMatrixMath.RotateRows(b, i, j, c, s);
            InternalMatrixMath.RotateColumns(q, i, j, c, C.Conjugate(s));
        }

        /// <summary>Applies a right rotation to both matrices and updates Z.</summary>
        /// <param name="a">First work matrix.</param>
        /// <param name="b">Second work matrix.</param>
        /// <param name="z">Right accumulator.</param>
        /// <param name="i">First column.</param>
        /// <param name="j">Second column.</param>
        /// <param name="c">Real cosine.</param>
        /// <param name="s">Complex sine acting on ordered columns i,j.</param>
        private static void RightPair(C[,] a, C[,] b, C[,] z, int i, int j, double c, C s)
        {
            InternalMatrixMath.RotateColumns(a, i, j, c, s);
            InternalMatrixMath.RotateColumns(b, i, j, c, s);
            InternalMatrixMath.RotateColumns(z, i, j, c, s);
        }

    }
}
