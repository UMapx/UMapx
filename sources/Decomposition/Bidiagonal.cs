using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides real and complex bidiagonal reduction.</summary>
    public static class Bidiagonal
    {
        /// <summary>Computes A = U B V^T by two-sided Householder reduction.</summary>
        /// <param name="matrix">Finite nonempty rectangular matrix.</param>
        /// <returns>Full square U and V and upper bidiagonal B with the same dimensions as the input.</returns>
        public static (float[,] U, float[,] B, float[,] V) Decompose(float[,] matrix)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix));
            return (InternalMatrixMath.Real(d.U), InternalMatrixMath.Real(d.B), InternalMatrixMath.Real(d.V));
        }

        /// <summary>Computes A = U B V^H by two-sided Householder reduction.</summary>
        /// <param name="matrix">Finite nonempty rectangular matrix.</param>
        /// <returns>Full square U and V and upper bidiagonal B with the same dimensions as the input.</returns>
        public static (Complex32[,] U, Complex32[,] B, Complex32[,] V) Decompose(Complex32[,] matrix)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix));
            return (InternalMatrixMath.Single(d.U), InternalMatrixMath.Single(d.B), InternalMatrixMath.Single(d.V));
        }

        /// <summary>Alternates left and right reflections to obtain upper bidiagonal form.</summary>
        /// <param name="a">Private rectangular input buffer, overwritten by B.</param>
        /// <returns>Accumulated unitary factors and the bidiagonal buffer.</returns>
        internal static (C[,] U, C[,] B, C[,] V) Factor(C[,] a)
        {
            int m = a.GetLength(0), n = a.GetLength(1);
            var u = InternalMatrixMath.Eye(m);
            var v = InternalMatrixMath.Eye(n);
            for (int k = 0; k < Math.Min(m, n); k++)
            {
                var left = InternalMatrixMath.Column(a, k, k);
                left = InternalMatrixMath.HouseholderVector(left);
                InternalMatrixMath.ReflectLeft(a, left, k, k);
                InternalMatrixMath.ReflectRight(u, left, k, 0);
                for (int i = k + 1; i < m; i++) a[i, k] = 0;
                if (k + 1 >= n) continue;
                var right = InternalMatrixMath.ConjugateRow(a, k, k + 1);
                right = InternalMatrixMath.HouseholderVector(right);
                InternalMatrixMath.ReflectRight(a, right, k + 1, k);
                InternalMatrixMath.ReflectRight(v, right, k + 1, 0);
                for (int j = k + 2; j < n; j++) a[k, j] = 0;
            }
            return (u, a, v);
        }
    }
}
