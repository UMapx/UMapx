using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides real and complex bidiagonal reduction</summary>
    public static class Bidiagonal
    {
        /// <summary>Computes A = U B V^T by two-sided Householder reduction</summary>
        /// <param name="matrix">Finite nonempty rectangular matrix.</param>
        /// <returns>Full square U and V and upper bidiagonal B with the same dimensions as the input.</returns>
        public static (float[,] U, float[,] B, float[,] V) Decompose(float[,] matrix)
        {
            var d = Factor(MatrixMath.Copy(matrix));
            return (MatrixMath.Real(d.U), MatrixMath.Real(d.B), MatrixMath.Real(d.V));
        }

        /// <summary>Computes A = U B V^H by two-sided Householder reduction</summary>
        /// <param name="matrix">Finite nonempty rectangular matrix.</param>
        /// <returns>Full square U and V and upper bidiagonal B with the same dimensions as the input.</returns>
        public static (Complex32[,] U, Complex32[,] B, Complex32[,] V) Decompose(Complex32[,] matrix)
        {
            var d = Factor(MatrixMath.Copy(matrix));
            return (MatrixMath.Single(d.U), MatrixMath.Single(d.B), MatrixMath.Single(d.V));
        }

        /// <summary>Alternates left and right reflections to obtain upper bidiagonal form</summary>
        /// <param name="a">Private rectangular input buffer, overwritten by B.</param>
        /// <returns>Accumulated unitary factors and the bidiagonal buffer.</returns>
        internal static (C[,] U, C[,] B, C[,] V) Factor(C[,] a)
        {
            int m = a.GetLength(0), n = a.GetLength(1);
            var u = MatrixMath.Eye(m);
            var v = MatrixMath.Eye(n);
            for (int k = 0; k < Math.Min(m, n); k++)
            {
                var left = new C[m - k];
                for (int i = k; i < m; i++) left[i - k] = a[i, k];
                left = Householder.Vector(left);
                Householder.ApplyLeft(a, left, k, k);
                Householder.ApplyRight(u, left, k, 0);
                for (int i = k + 1; i < m; i++) a[i, k] = 0;
                if (k + 1 >= n) continue;
                var right = new C[n - k - 1];
                for (int j = k + 1; j < n; j++) right[j - k - 1] = C.Conjugate(a[k, j]);
                right = Householder.Vector(right);
                Householder.ApplyRight(a, right, k + 1, k);
                Householder.ApplyRight(v, right, k + 1, 0);
                for (int j = k + 2; j < n; j++) a[k, j] = 0;
            }
            return (u, a, v);
        }
    }
}
