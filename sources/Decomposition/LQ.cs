using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides real and complex LQ decomposition.</summary>
    public static class LQ
    {
        /// <summary>Computes the economy-size LQ factorization of a rectangular matrix.</summary>
        /// <param name="matrix">Finite nonempty input matrix, not modified.</param>
        /// <returns>Factors L, Q; their product equals the input. Q has orthonormal rows.</returns>
        public static (float[,] L, float[,] Q) Decompose(float[,] matrix)
        {
            MatrixMath.CheckShape(matrix);
            var d = QR.Decompose(matrix.Transpose());
            return (d.R.Transpose(), d.Q.Transpose());
        }

        /// <summary>Computes the economy-size LQ factorization of a rectangular matrix.</summary>
        /// <param name="matrix">Finite nonempty input matrix, not modified.</param>
        /// <returns>Factors L, Q; their product equals the input. Q has orthonormal rows.</returns>
        public static (Complex32[,] L, Complex32[,] Q) Decompose(Complex32[,] matrix)
        {
            MatrixMath.CheckShape(matrix);
            var d = QR.Decompose(matrix.Hermitian());
            return (d.R.Hermitian(), d.Q.Hermitian());
        }


    }
}
