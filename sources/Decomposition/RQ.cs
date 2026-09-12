using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides real and complex RQ decomposition</summary>
    public static class RQ
    {
        /// <summary>Computes the economy-size RQ factorization of a rectangular matrix</summary>
        /// <param name="matrix">Finite nonempty input matrix, not modified.</param>
        /// <returns>Factors R, Q; their product equals the input. Q has orthonormal rows.</returns>
        public static (float[,] R, float[,] Q) Decompose(float[,] matrix)
        {
            MatrixMath.CheckShape(matrix);
            var d = QR.Decompose(matrix.Flip(Direction.Vertical).Transpose());
            return (d.R.Transpose().Flip(Direction.Both), d.Q.Transpose().Flip(Direction.Vertical));
        }

        /// <summary>Computes the economy-size RQ factorization of a rectangular matrix</summary>
        /// <param name="matrix">Finite nonempty input matrix, not modified.</param>
        /// <returns>Factors R, Q; their product equals the input. Q has orthonormal rows.</returns>
        public static (Complex32[,] R, Complex32[,] Q) Decompose(Complex32[,] matrix)
        {
            MatrixMath.CheckShape(matrix);
            var d = QR.Decompose(matrix.Flip(Direction.Vertical).Hermitian());
            return (d.R.Hermitian().Flip(Direction.Both), d.Q.Hermitian().Flip(Direction.Vertical));
        }


    }
}
