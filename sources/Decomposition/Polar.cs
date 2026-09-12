using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides right polar decomposition for real and complex rectangular matrices</summary>
    public static class Polar
    {
        /// <summary>Computes the right polar decomposition A = U P</summary>
        /// <param name="matrix">Finite nonempty rectangular matrix.</param>
        /// <param name="iterations">Positive SVD iteration limit.</param>
        /// <returns>The partial isometry U and positive semidefinite symmetric P. Square full-rank U is orthogonal.</returns>
        public static (float[,] U, float[,] P) Decompose(float[,] matrix, int iterations = 10)
        {
            var d = SVD.Decompose(matrix, iterations);
            var right = d.V.Transpose();
            return (d.U.Dot(right), d.V.Dot(d.S).Dot(right));
        }

        /// <summary>Computes the right polar decomposition A = U P</summary>
        /// <param name="matrix">Finite nonempty rectangular matrix.</param>
        /// <param name="iterations">Positive SVD iteration limit.</param>
        /// <returns>The partial isometry U and positive semidefinite Hermitian P. Square full-rank U is unitary.</returns>
        public static (Complex32[,] U, Complex32[,] P) Decompose(Complex32[,] matrix, int iterations = 50)
        {
            var d = SVD.Decompose(matrix, iterations);
            var right = d.V.Hermitian();
            return (d.U.Dot(right), d.V.Dot(d.S).Dot(right));
        }


    }
}
