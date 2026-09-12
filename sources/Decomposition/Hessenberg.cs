using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides orthogonal and unitary Hessenberg reductions.</summary>
    public static class Hessenberg
    {
        /// <summary>Computes A = P H P^T.</summary>
        /// <param name="matrix">Finite nonempty square matrix, not modified.</param>
        /// <returns>Orthogonal P and upper Hessenberg H.</returns>
        public static (float[,] P, float[,] H) Decompose(float[,] matrix)
        {
            var d = Factor(InternalMatrixMath.CopyReal(matrix, true));
            return (InternalMatrixMath.Real(d.P), InternalMatrixMath.Real(d.H));
        }

        /// <summary>Computes A = P H P^H.</summary>
        /// <param name="matrix">Finite nonempty complex square matrix, not modified.</param>
        /// <returns>Unitary P and upper Hessenberg H.</returns>
        public static (Complex32[,] P, Complex32[,] H) Decompose(Complex32[,] matrix)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix, true));
            return (InternalMatrixMath.Single(d.P), InternalMatrixMath.Single(d.H));
        }

        /// <summary>Applies two-sided Householder similarities in double precision.</summary>
        /// <param name="a">Private square work buffer, overwritten by Hessenberg form.</param>
        /// <returns>The full unitary accumulator and reduced buffer.</returns>
        internal static (C[,] P, C[,] H) Factor(C[,] a)
        {
            int n = a.GetLength(0);
            var p = InternalMatrixMath.Eye(n);
            for (int k = 0; k < n - 2; k++)
            {
                var v = InternalMatrixMath.Column(a, k, k + 1);
                v = InternalMatrixMath.HouseholderVector(v);
                InternalMatrixMath.ReflectLeft(a, v, k + 1, k);
                InternalMatrixMath.ReflectRight(a, v, k + 1, 0);
                InternalMatrixMath.ReflectRight(p, v, k + 1, 0);
                for (int i = k + 2; i < n; i++) a[i, k] = 0;
            }
            return (p, a);
        }

        /// <summary>Applies two-sided Householder similarities in double precision.</summary>
        /// <param name="a">Private square work buffer, overwritten by Hessenberg form.</param>
        /// <returns>The full orthogonal accumulator and reduced buffer.</returns>
        internal static (double[][] P, double[][] H) Factor(double[][] a)
        {
            int n = a.Length;
            var scratch = new double[n];
            var p = InternalMatrixMath.EyeJagged(n);
            for (int k = 0; k < n - 2; k++)
            {
                var v = InternalMatrixMath.Column(a, k, k + 1);
                v = InternalMatrixMath.HouseholderVector(v);
                InternalMatrixMath.ReflectLeft(a, v, k + 1, k, scratch);
                InternalMatrixMath.ReflectRight(a, v, k + 1, 0);
                InternalMatrixMath.ReflectRight(p, v, k + 1, 0);
                for (int i = k + 2; i < n; i++) a[i][k] = 0;
            }
            return (p, a);
        }
    }
}
