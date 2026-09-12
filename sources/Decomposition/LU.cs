using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides LU decomposition with partial row pivoting.</summary>
    public static class LU
    {
        /// <summary>Computes A[P,:] = L U with partial row pivoting.</summary>
        /// <param name="matrix">Finite nonempty square matrix, not modified.</param>
        /// <returns>Unit lower triangular L, upper triangular U, and the row permutation P.</returns>
        public static (float[,] L, float[,] U, int[] P) Decompose(float[,] matrix)
        {
            var d = Factor(InternalMatrixMath.CopyReal(matrix, true));
            return (InternalMatrixMath.Real(d.L), InternalMatrixMath.Real(d.U), d.P);
        }

        /// <summary>Computes A[P,:] = L U with complex partial row pivoting.</summary>
        /// <param name="matrix">Finite nonempty square matrix, not modified.</param>
        /// <returns>Unit lower triangular L, upper triangular U, and the row permutation P.</returns>
        public static (Complex32[,] L, Complex32[,] U, int[] P) Decompose(Complex32[,] matrix)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix, true));
            return (InternalMatrixMath.Single(d.L), InternalMatrixMath.Single(d.U), d.P);
        }

        /// <summary>Builds the row permutation matrix from a pivot vector.</summary>
        /// <param name="permutation">A permutation of indices zero through n-1.</param>
        /// <returns>P such that P*A selects rows A[permutation[i],:].</returns>
        public static float[,] PermutationMatrix(int[] permutation)
        {
            if (permutation == null) throw new ArgumentNullException(nameof(permutation));
            int n = permutation.Length;
            var p = new float[n, n];
            var used = new bool[n];
            for (int i = 0; i < n; i++)
            {
                int j = permutation[i];
                if (j < 0 || j >= n || used[j]) throw new ArgumentException("Invalid row permutation.", nameof(permutation));
                used[j] = true; p[i, j] = 1;
            }
            return p;
        }

        /// <summary>Performs Gaussian elimination with magnitude-based row pivoting.</summary>
        /// <param name="a">Private square buffer overwritten by U.</param>
        /// <returns>L, U and row indices, including valid factors of singular matrices.</returns>
        internal static (C[,] L, C[,] U, int[] P) Factor(C[,] a)
        {
            int n = a.GetLength(0);
            var l = InternalMatrixMath.Eye(n);
            var p = new int[n];
            for (int i = 0; i < n; i++) p[i] = i;
            for (int k = 0; k < n; k++)
            {
                int pivot = k;
                for (int i = k + 1; i < n; i++)
                    if (C.Abs(a[i, k]) > C.Abs(a[pivot, k])) pivot = i;
                if (pivot != k)
                {
                    int t = p[k]; p[k] = p[pivot]; p[pivot] = t;
                    for (int j = 0; j < n; j++) { C z = a[k, j]; a[k, j] = a[pivot, j]; a[pivot, j] = z; }
                    for (int j = 0; j < k; j++) { C z = l[k, j]; l[k, j] = l[pivot, j]; l[pivot, j] = z; }
                }
                if (a[k, k] == C.Zero) continue;
                for (int i = k + 1; i < n; i++)
                {
                    l[i, k] = a[i, k] / a[k, k];
                    a[i, k] = 0;
                    for (int j = k + 1; j < n; j++) a[i, j] -= l[i, k] * a[k, j];
                }
            }
            return (l, a, p);
        }


        /// <summary>Performs Gaussian elimination with magnitude-based row pivoting.</summary>
        /// <param name="a">Private square buffer overwritten by U.</param>
        /// <returns>L, U and row indices, including valid factors of singular matrices.</returns>
        internal static (double[][] L, double[][] U, int[] P) Factor(double[][] a)
        {
            int n = a.Length;
            var l = InternalMatrixMath.EyeJagged(n);
            var p = new int[n];
            for (int i = 0; i < n; i++) p[i] = i;
            for (int k = 0; k < n; k++)
            {
                int pivot = k;
                for (int i = k + 1; i < n; i++)
                    if (Math.Abs(a[i][k]) > Math.Abs(a[pivot][k])) pivot = i;
                if (pivot != k)
                {
                    int t = p[k]; p[k] = p[pivot]; p[pivot] = t;
                    for (int j = 0; j < n; j++) { double z = a[k][j]; a[k][j] = a[pivot][j]; a[pivot][j] = z; }
                    for (int j = 0; j < k; j++) { double z = l[k][j]; l[k][j] = l[pivot][j]; l[pivot][j] = z; }
                }
                if (a[k][k] == 0) continue;
                for (int i = k + 1; i < n; i++)
                {
                    l[i][k] = a[i][k] / a[k][k];
                    a[i][k] = 0;
                    for (int j = k + 1; j < n; j++) a[i][j] -= l[i][k] * a[k][j];
                }
            }
            return (l, a, p);
        }
    }
}
