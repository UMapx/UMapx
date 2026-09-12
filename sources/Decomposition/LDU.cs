using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides pivoted LDU decomposition.</summary>
    public static class LDU
    {
        /// <summary>Computes A[P,:] = L diag(D) U with unit triangular factors.</summary>
        /// <param name="matrix">Finite nonempty square matrix with nonzero LU pivots.</param>
        /// <returns>L, diagonal D, U, and the row permutation P.</returns>
        public static (float[,] L, float[] D, float[,] U, int[] P) Decompose(float[,] matrix)
        {
            var lu = LU.Decompose(matrix);
            int n = matrix.GetLength(0);
            var d = new float[n];
            for (int i = 0; i < n; i++)
            {
                d[i] = lu.U[i, i];
                if (d[i] == 0) throw new InvalidOperationException("A zero pivot prevents unit-diagonal LDU factorization.");
                for (int j = i; j < n; j++) lu.U[i, j] /= d[i];
            }
            return (lu.L, d, lu.U, lu.P);
        }

        /// <summary>Computes A[P,:] = L diag(D) U with unit triangular factors.</summary>
        /// <param name="matrix">Finite nonempty square matrix with nonzero LU pivots.</param>
        /// <returns>L, diagonal D, U, and the row permutation P.</returns>
        public static (Complex32[,] L, Complex32[] D, Complex32[,] U, int[] P) Decompose(Complex32[,] matrix)
        {
            var lu = LU.Decompose(matrix);
            int n = matrix.GetLength(0);
            var d = new Complex32[n];
            for (int i = 0; i < n; i++)
            {
                d[i] = lu.U[i, i];
                if (d[i] == 0) throw new InvalidOperationException("A zero pivot prevents unit-diagonal LDU factorization.");
                for (int j = i; j < n; j++) lu.U[i, j] /= d[i];
            }
            return (lu.L, d, lu.U, lu.P);
        }


    }
}
