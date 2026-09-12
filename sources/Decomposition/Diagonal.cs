using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides diagonal extraction factorization.</summary>
    public static class Diagonal
    {
        /// <summary>Computes A = B diag(D) by extracting and dividing by the diagonal.</summary>
        /// <param name="matrix">Finite nonempty square matrix with nonzero diagonal entries.</param>
        /// <returns>A column-normalized B and the original diagonal D.</returns>
        public static (float[,] B, float[] D) Decompose(float[,] matrix)
        {
            InternalRealMatrixMath.Validate(matrix, true);
            int n = matrix.GetLength(0);
            var d = new float[n];
            var b = (float[,])matrix.Clone();
            for (int j = 0; j < n; j++)
            {
                d[j] = matrix[j, j];
                if (d[j] == 0) throw new ArgumentException("Diagonal entries must be nonzero.", nameof(matrix));
                for (int i = 0; i < n; i++) b[i, j] /= d[j];
            }
            return (b, d);
        }

        /// <summary>Computes A = B diag(D) by extracting and dividing by the diagonal.</summary>
        /// <param name="matrix">Finite nonempty square matrix with nonzero diagonal entries.</param>
        /// <returns>A column-normalized B and the original diagonal D.</returns>
        public static (Complex32[,] B, Complex32[] D) Decompose(Complex32[,] matrix)
        {
            InternalMatrixMath.Copy(matrix, true);
            int n = matrix.GetLength(0);
            var d = new Complex32[n];
            var b = (Complex32[,])matrix.Clone();
            for (int j = 0; j < n; j++)
            {
                d[j] = matrix[j, j];
                if (d[j] == 0) throw new ArgumentException("Diagonal entries must be nonzero.", nameof(matrix));
                for (int i = 0; i < n; i++) b[i, j] /= d[j];
            }
            return (b, d);
        }


    }
}
