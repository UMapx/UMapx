using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides real and complex QL decomposition</summary>
    public static class QL
    {
        /// <summary>Computes the economy-size QL factorization of a rectangular matrix</summary>
        /// <param name="matrix">Finite nonempty input matrix, not modified.</param>
        /// <returns>Factors Q, L; their product equals the input. Q has orthonormal columns.</returns>
        public static (float[,] Q, float[,] L) Decompose(float[,] matrix)
        {
            InternalMatrixMath.CheckShape(matrix);
            var d = QR.Decompose(matrix.Flip(Direction.Horizontal));
            return (d.Q.Flip(Direction.Horizontal), d.R.Flip(Direction.Both));
        }

        /// <summary>Computes the economy-size QL factorization of a rectangular matrix</summary>
        /// <param name="matrix">Finite nonempty input matrix, not modified.</param>
        /// <returns>Factors Q, L; their product equals the input. Q has orthonormal columns.</returns>
        public static (Complex32[,] Q, Complex32[,] L) Decompose(Complex32[,] matrix)
        {
            InternalMatrixMath.CheckShape(matrix);
            var d = QR.Decompose(matrix.Flip(Direction.Horizontal));
            return (d.Q.Flip(Direction.Horizontal), d.R.Flip(Direction.Both));
        }


    }
}
