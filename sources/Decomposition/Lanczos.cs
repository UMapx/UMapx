using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides symmetric and Hermitian Lanczos reductions with orthogonal restarts.</summary>
    public static class Lanczos
    {
        /// <summary>Computes a symmetric Lanczos tridiagonalization.</summary>
        /// <param name="matrix">Finite nonempty symmetric square matrix.</param>
        /// <param name="full">Whether to always perform a second reorthogonalization pass; otherwise it is applied when cancellation requires it.</param>
        /// <returns>Q and tridiagonal T satisfying A = Q T Q^T.</returns>
        public static (float[,] Q, float[,] T) Decompose(float[,] matrix, bool full = false)
        {
            var a = InternalMatrixMath.CopyReal(matrix, true);
            InternalMatrixMath.RequireSymmetric(a);
            var d = Arnoldi.Factor(a, full);
            int n = a.GetLength(0);
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++) d.H[i][j] = j == i + 1 ? d.H[j][i] : 0;
            }
            return (InternalMatrixMath.Real(d.Q), InternalMatrixMath.Real(d.H));
        }

        /// <summary>Computes a Hermitian Lanczos tridiagonalization.</summary>
        /// <param name="matrix">Finite nonempty Hermitian square matrix.</param>
        /// <param name="full">Whether to always perform a second reorthogonalization pass; otherwise it is applied when cancellation requires it.</param>
        /// <returns>Q and tridiagonal T satisfying A = Q T Q^H.</returns>
        public static (Complex32[,] Q, Complex32[,] T) Decompose(Complex32[,] matrix, bool full = false)
        {
            var a = InternalMatrixMath.Copy(matrix, true);
            InternalMatrixMath.RequireHermitian(a);
            var d = Arnoldi.Factor(a, full);
            int n = a.GetLength(0);
            for (int i = 0; i < n; i++)
            {
                d.H[i, i] = d.H[i, i].Real;
                for (int j = i + 1; j < n; j++) d.H[i, j] = j == i + 1 ? C.Conjugate(d.H[j, i]) : C.Zero;
            }
            return (InternalMatrixMath.Single(d.Q), InternalMatrixMath.Single(d.H));
        }


    }
}
