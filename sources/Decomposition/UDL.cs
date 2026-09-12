using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides unpivoted UDL factorization for symmetric and Hermitian matrices.</summary>
    public static class UDL
    {
        /// <summary>Computes A = U diag(D) U^T without diagonal pivoting.</summary>
        /// <param name="matrix">Finite nonempty symmetric matrix with nonzero elimination pivots.</param>
        /// <returns>Unit upper triangular U and real diagonal D. Indefinite inputs are supported when no pivot vanishes.</returns>
        public static (float[,] U, float[] D) Decompose(float[,] matrix)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix, true));
            return (InternalMatrixMath.Real(d.F), d.D);
        }

        /// <summary>Constructs the conjugate-transposed factor from an existing factor.</summary>
        /// <param name="factor">Square triangular factor from Decompose.</param>
        /// <returns>The lower factor.</returns>
        public static float[,] LowerFactor(float[,] factor) => InternalMatrixMath.Real(InternalMatrixMath.Adjoint(InternalMatrixMath.Copy(factor, true)));

        /// <summary>Computes A = U diag(D) U^H without diagonal pivoting.</summary>
        /// <param name="matrix">Finite nonempty Hermitian matrix with nonzero elimination pivots.</param>
        /// <returns>Unit upper triangular U and real diagonal D. Indefinite inputs are supported when no pivot vanishes.</returns>
        public static (Complex32[,] U, float[] D) Decompose(Complex32[,] matrix)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix, true));
            return (InternalMatrixMath.Single(d.F), d.D);
        }

        /// <summary>Constructs the conjugate-transposed factor from an existing factor.</summary>
        /// <param name="factor">Square triangular factor from Decompose.</param>
        /// <returns>The lower factor.</returns>
        public static Complex32[,] LowerFactor(Complex32[,] factor) => InternalMatrixMath.Single(InternalMatrixMath.Adjoint(InternalMatrixMath.Copy(factor, true)));

        /// <summary>Performs Hermitian diagonal elimination in double precision without pivoting.</summary>
        /// <param name="a">Private Hermitian square buffer.</param>
        /// <returns>A unit triangular factor and real diagonal; a zero pivot is rejected.</returns>
        private static (C[,] F, float[] D) Factor(C[,] a)
        {
            InternalMatrixMath.RequireHermitian(a);
            int n = a.GetLength(0);
            var f = InternalMatrixMath.Eye(n);
            var d = new double[n];
            for (int step = 0; step < n; step++)
            {
                int j = n - 1 - step;
                double pivot = a[j, j].Real;
                for (int prev = 0; prev < step; prev++)
                {
                    int k = n - 1 - prev;
                    double magnitude = C.Abs(f[j, k]);
                    pivot -= magnitude * magnitude * d[k];
                }
                if (pivot == 0) throw new InvalidOperationException("A zero pivot requires a pivoted Hermitian factorization.");
                d[j] = pivot;
                for (int next = step + 1; next < n; next++)
                {
                    int i = n - 1 - next;
                    C sum = a[i, j];
                    for (int prev = 0; prev < step; prev++)
                    {
                        int k = n - 1 - prev;
                        sum -= f[i, k] * d[k] * C.Conjugate(f[j, k]);
                    }
                    f[i, j] = sum / pivot;
                }
            }
            var diagonal = new float[n];
            for (int i = 0; i < n; i++) diagonal[i] = (float)d[i];
            return (f, diagonal);
        }
    }
}
