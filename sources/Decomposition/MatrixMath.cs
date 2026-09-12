using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides shared double-precision work-buffer operations for matrix decompositions.</summary>
    internal static class MatrixMath
    {
        internal const double Roundoff = 2.2204460492503131e-16;
        internal const double SingleRoundoff = 1.1920928955078125e-7;

        /// <summary>Copies a finite, nonempty real matrix into complex double-precision storage.</summary>
        /// <param name="a">Input matrix, which is not modified.</param>
        /// <param name="square">Whether equal dimensions are required.</param>
        /// <returns>An independent work buffer.</returns>
        internal static C[,] Copy(float[,] a, bool square = false)
        {
            CheckShape(a, square);
            var b = new C[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < b.GetLength(0); i++)
                for (int j = 0; j < b.GetLength(1); j++)
                {
                    if (float.IsNaN(a[i, j]) || float.IsInfinity(a[i, j]))
                        throw new ArgumentException("The matrix must contain finite values.", nameof(a));
                    b[i, j] = a[i, j];
                }
            return b;
        }

        /// <summary>Copies a finite, nonempty complex matrix into double-precision storage.</summary>
        /// <param name="a">Input matrix, which is not modified.</param>
        /// <param name="square">Whether equal dimensions are required.</param>
        /// <returns>An independent work buffer.</returns>
        internal static C[,] Copy(Complex32[,] a, bool square = false)
        {
            CheckShape(a, square);
            var b = new C[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < b.GetLength(0); i++)
                for (int j = 0; j < b.GetLength(1); j++)
                {
                    var z = a[i, j];
                    if (float.IsNaN(z.Real) || float.IsInfinity(z.Real) || float.IsNaN(z.Imag) || float.IsInfinity(z.Imag))
                        throw new ArgumentException("The matrix must contain finite values.", nameof(a));
                    b[i, j] = new C(z.Real, z.Imag);
                }
            return b;
        }

        /// <summary>Validates matrix dimensions before allocating numerical work buffers.</summary>
        /// <param name="a">A two-dimensional array.</param>
        /// <param name="square">Whether a square matrix is required.</param>
        internal static void CheckShape(Array a, bool square = false)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (a.Rank != 2 || a.GetLength(0) == 0 || a.GetLength(1) == 0)
                throw new ArgumentException("The matrix must be nonempty and two-dimensional.", nameof(a));
            if (square && a.GetLength(0) != a.GetLength(1))
                throw new ArgumentException("The matrix must be square.", nameof(a));
        }

        /// <summary>Narrows a work matrix to complex single precision.</summary>
        /// <param name="a">Double-precision values.</param>
        /// <returns>A newly allocated complex matrix.</returns>
        internal static Complex32[,] Single(C[,] a)
        {
            var b = new Complex32[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < b.GetLength(0); i++)
                for (int j = 0; j < b.GetLength(1); j++) b[i, j] = new Complex32((float)a[i, j].Real, (float)a[i, j].Imaginary);
            return b;
        }

        /// <summary>Narrows a real-valued work matrix to single precision.</summary>
        /// <param name="a">Work buffer whose imaginary components are zero.</param>
        /// <returns>A newly allocated real matrix.</returns>
        internal static float[,] Real(C[,] a)
        {
            var b = new float[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < b.GetLength(0); i++)
                for (int j = 0; j < b.GetLength(1); j++) b[i, j] = (float)a[i, j].Real;
            return b;
        }

        /// <summary>Creates an identity matrix for accumulating unitary transformations.</summary>
        /// <param name="n">Nonnegative order.</param>
        /// <returns>The identity of order n.</returns>
        internal static C[,] Eye(int n)
        {
            var a = new C[n, n];
            for (int i = 0; i < n; i++) a[i, i] = C.One;
            return a;
        }

        /// <summary>Computes conjugate transposition, including for real-valued work buffers.</summary>
        /// <param name="a">Input matrix.</param>
        /// <returns>The conjugate transpose.</returns>
        internal static C[,] Adjoint(C[,] a)
        {
            var b = new C[a.GetLength(1), a.GetLength(0)];
            for (int i = 0; i < a.GetLength(0); i++)
                for (int j = 0; j < a.GetLength(1); j++) b[j, i] = C.Conjugate(a[i, j]);
            return b;
        }

        /// <summary>Multiplies compatible work matrices with complex double accumulation.</summary>
        /// <param name="a">Left matrix.</param>
        /// <param name="b">Right matrix.</param>
        /// <returns>The matrix product.</returns>
        internal static C[,] Multiply(C[,] a, C[,] b)
        {
            if (a.GetLength(1) != b.GetLength(0)) throw new ArgumentException("Incompatible matrix dimensions.");
            var c = new C[a.GetLength(0), b.GetLength(1)];
            for (int i = 0; i < c.GetLength(0); i++)
                for (int k = 0; k < a.GetLength(1); k++)
                    for (int j = 0; j < c.GetLength(1); j++) c[i, j] += a[i, k] * b[k, j];
            return c;
        }

        /// <summary>Computes a scaled Euclidean norm without squaring large or tiny entries directly.</summary>
        /// <param name="v">Vector with finite entries.</param>
        /// <returns>The nonnegative Euclidean norm.</returns>
        internal static double Norm(C[] v)
        {
            double scale = 0, sum = 1;
            foreach (C value in v)
            {
                double x = C.Abs(value);
                if (x == 0) continue;
                if (scale < x) { double r = scale / x; sum = 1 + sum * r * r; scale = x; }
                else { double r = x / scale; sum += r * r; }
            }
            return scale == 0 ? 0 : scale * Math.Sqrt(sum);
        }

        /// <summary>Returns the largest entry magnitude, or zero for a zero matrix.</summary>
        /// <param name="a">Finite work matrix.</param>
        /// <returns>A nonnegative scale.</returns>
        internal static double Max(C[,] a)
        {
            double scale = 0;
            foreach (C z in a) scale = Math.Max(scale, C.Abs(z));
            return scale;
        }

        /// <summary>Checks the Hermitian condition with a relative single-precision tolerance.</summary>
        /// <param name="a">Square input matrix.</param>
        internal static void RequireHermitian(C[,] a)
        {
            if (!IsHermitian(a)) throw new ArgumentException("The matrix must be Hermitian (symmetric for real inputs).");
        }

        /// <summary>Detects Hermitian structure to the resolution of single-precision input.</summary>
        /// <param name="a">Finite square work matrix.</param>
        /// <returns>True when conjugate symmetry holds within eight single-precision rounding units of the matrix scale.</returns>
        internal static bool IsHermitian(C[,] a)
        {
            int n = a.GetLength(0);
            if (n != a.GetLength(1)) return false;
            double tolerance = 8 * SingleRoundoff * Max(a);
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= i; j++)
                    if (C.Abs(a[i, j] - C.Conjugate(a[j, i])) > tolerance)
                        return false;
            return true;
        }

        /// <summary>Copies a rectangular leading block or a block starting at specified offsets.</summary>
        /// <param name="a">Source matrix.</param>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        /// <param name="row">First source row.</param>
        /// <param name="column">First source column.</param>
        /// <returns>An independent block.</returns>
        internal static C[,] Block(C[,] a, int rows, int columns, int row = 0, int column = 0)
        {
            var b = new C[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++) b[i, j] = a[i + row, j + column];
            return b;
        }

        /// <summary>Completes an orthonormal basis after exact or numerical breakdown.</summary>
        /// <param name="q">Matrix containing orthonormal columns before column k.</param>
        /// <param name="k">Column to complete.</param>
        /// <returns>A unit vector orthogonal to the preceding columns.</returns>
        internal static C[] Complete(C[,] q, int k)
        {
            int n = q.GetLength(0);
            var best = new C[n];
            double largest = -1;
            for (int seed = 0; seed < n; seed++)
            {
                var v = new C[n]; v[seed] = 1;
                Orthogonalize(v, q, k);
                double norm = Norm(v);
                if (norm > largest) { largest = norm; best = v; }
            }
            if (largest <= Roundoff) throw new InvalidOperationException("Unable to complete an orthonormal basis.");
            for (int i = 0; i < n; i++) best[i] /= largest;
            return best;
        }

        /// <summary>Applies two passes of modified Gram-Schmidt using the Hermitian inner product.</summary>
        /// <param name="v">Vector modified in place.</param>
        /// <param name="q">Previously computed orthonormal columns.</param>
        /// <param name="columns">Number of columns to remove.</param>
        internal static void Orthogonalize(C[] v, C[,] q, int columns)
        {
            for (int pass = 0; pass < 2; pass++)
                for (int j = 0; j < columns; j++)
                {
                    C dot = 0;
                    for (int i = 0; i < v.Length; i++) dot += C.Conjugate(q[i, j]) * v[i];
                    for (int i = 0; i < v.Length; i++) v[i] -= q[i, j] * dot;
                }
        }
    }
}
