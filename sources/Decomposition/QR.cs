using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides Householder QR factorization and access to its reflection vectors</summary>
    public static class QR
    {
        /// <summary>Computes A = Q R using an economy-size orthogonal factor</summary>
        /// <param name="matrix">Finite nonempty m by n matrix; not modified.</param>
        /// <returns>Q of size m by min(m,n) and R of size min(m,n) by n.</returns>
        public static (float[,] Q, float[,] R) Decompose(float[,] matrix)
        {
            var d = Factor(MatrixMath.Copy(matrix), full: false);
            int k = Math.Min(matrix.GetLength(0), matrix.GetLength(1));
            return (MatrixMath.Real(MatrixMath.Block(d.Q, matrix.GetLength(0), k)),
                    MatrixMath.Real(MatrixMath.Block(d.R, k, matrix.GetLength(1))));
        }

        /// <summary>Computes A = Q R with Q^H Q = I using complex Householder reflections</summary>
        /// <param name="matrix">Finite nonempty m by n matrix; not modified.</param>
        /// <returns>Q of size m by min(m,n) and R of size min(m,n) by n.</returns>
        public static (Complex32[,] Q, Complex32[,] R) Decompose(Complex32[,] matrix)
        {
            var d = Factor(MatrixMath.Copy(matrix), full: false);
            int k = Math.Min(matrix.GetLength(0), matrix.GetLength(1));
            return (MatrixMath.Single(MatrixMath.Block(d.Q, matrix.GetLength(0), k)),
                    MatrixMath.Single(MatrixMath.Block(d.R, k, matrix.GetLength(1))));
        }

        /// <summary>Computes normalized reflection vectors without constructing Q</summary>
        /// <param name="matrix">Finite input matrix; a separate reduction is performed.</param>
        /// <returns>An m by min(m,n) matrix of vectors v defining H = I - 2 v v^T; zero columns denote identity.</returns>
        public static float[,] HouseholderVectors(float[,] matrix) => MatrixMath.Real(Factor(MatrixMath.Copy(matrix), false).H);

        /// <summary>Computes normalized complex reflection vectors without constructing Q</summary>
        /// <param name="matrix">Finite input matrix; a separate reduction is performed.</param>
        /// <returns>An m by min(m,n) matrix of vectors v defining H = I - 2 v v^H; zero columns denote identity.</returns>
        public static Complex32[,] HouseholderVectors(Complex32[,] matrix) => MatrixMath.Single(Factor(MatrixMath.Copy(matrix), false).H);

        /// <summary>Reduces a private work matrix to upper trapezoidal form</summary>
        /// <param name="a">Work buffer overwritten by R.</param>
        /// <param name="vectors">Whether to construct Q from the reflection vectors.</param>
        /// <param name="full">Whether Q should be square instead of economy size.</param>
        /// <returns>Q (or null), R, and normalized reflection vectors.</returns>
        internal static (C[,] Q, C[,] R, C[,] H) Factor(C[,] a, bool vectors = true, bool full = true)
        {
            int m = a.GetLength(0), n = a.GetLength(1), kmax = Math.Min(m, n);
            C[,] q = null;
            var h = new C[m, kmax];
            for (int k = 0; k < kmax; k++)
            {
                var v = MatrixMath.Column(a, k, k);
                v = MatrixMath.HouseholderVector(v);
                MatrixMath.ReflectLeft(a, v, k, k);
                MatrixMath.SetColumn(h, k, v, k);
                for (int i = k + 1; i < m; i++) a[i, k] = 0;
            }
            if (vectors)
            {
                int columns = full ? m : kmax;
                q = new C[m, columns];
                for (int i = 0; i < columns; i++) q[i, i] = 1;
                // Reverse application builds Q without allocating an m by m matrix for a tall economy QR.
                for (int k = kmax - 1; k >= 0; k--)
                {
                    var v = MatrixMath.Column(h, k, k);
                    MatrixMath.ReflectLeft(q, v, k, 0);
                }
            }
            return (q, a, h);
        }
    }
}
