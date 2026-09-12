using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides Cholesky factorization for symmetric and Hermitian positive definite matrices.</summary>
    public static class Cholesky
    {
        /// <summary>Computes the lower triangular factor in A = L L^T.</summary>
        /// <param name="matrix">Finite nonempty symmetric positive definite matrix.</param>
        /// <returns>L with a strictly positive diagonal.</returns>
        public static float[,] Decompose(float[,] matrix) => InternalRealMatrixMath.Real(Factor(InternalRealMatrixMath.Copy(matrix, true)));

        /// <summary>Computes the lower triangular factor in A = L L^H.</summary>
        /// <param name="matrix">Finite nonempty Hermitian positive definite matrix.</param>
        /// <returns>L with a strictly positive real diagonal.</returns>
        public static Complex32[,] Decompose(Complex32[,] matrix) => InternalMatrixMath.Single(Factor(InternalMatrixMath.Copy(matrix, true)));

        /// <summary>Constructs the upper factor from a previously computed lower factor.</summary>
        /// <param name="lower">Square lower Cholesky factor.</param>
        /// <returns>L^T.</returns>
        public static float[,] UpperFactor(float[,] lower) => InternalRealMatrixMath.Real(InternalRealMatrixMath.Transpose(InternalRealMatrixMath.Copy(lower, true)));

        /// <summary>Constructs the upper factor from a previously computed complex lower factor.</summary>
        /// <param name="lower">Square lower Cholesky factor.</param>
        /// <returns>L^H.</returns>
        public static Complex32[,] UpperFactor(Complex32[,] lower) => InternalMatrixMath.Single(InternalMatrixMath.Adjoint(InternalMatrixMath.Copy(lower, true)));

        /// <summary>Computes Cholesky factors with Hermitian inner products in double precision.</summary>
        /// <param name="a">Private Hermitian square input.</param>
        /// <returns>A lower triangular factor; nonpositive pivots cause an exception.</returns>
        internal static C[,] Factor(C[,] a)
        {
            InternalMatrixMath.RequireHermitian(a);
            int n = a.GetLength(0);
            var l = new C[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= i; j++)
                {
                    C sum = a[i, j];
                    for (int k = 0; k < j; k++) sum -= l[i, k] * C.Conjugate(l[j, k]);
                    if (i == j)
                    {
                        if (!(sum.Real > 0)) throw new ArgumentException("The matrix must be positive definite.");
                        l[i, i] = Math.Sqrt(sum.Real);
                    }
                    else l[i, j] = sum / l[j, j].Real;
                }
            return l;
        }

        /// <summary>Computes Cholesky factors with symmetric inner products in double precision.</summary>
        /// <param name="a">Private symmetric square input.</param>
        /// <returns>A lower triangular factor; nonpositive pivots cause an exception.</returns>
        internal static double[][] Factor(double[][] a)
        {
            InternalRealMatrixMath.RequireSymmetric(a);
            int n = a.Length;
            var l = InternalRealMatrixMath.Create(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= i; j++)
                {
                    double sum = a[i][j] - InternalRealMatrixMath.Dot(l[i], l[j], j);
                    if (i == j)
                    {
                        if (!(sum > 0)) throw new ArgumentException("The matrix must be positive definite.");
                        l[i][i] = Math.Sqrt(sum);
                    }
                    else l[i][j] = sum / l[j][j];
                }
            return l;
        }
    }
}
