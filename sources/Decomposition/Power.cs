using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides real and complex power iteration.</summary>
    public static class Power
    {
        /// <summary>Approximates a dominant right eigenpair by normalized power iteration.</summary>
        /// <param name="matrix">Finite nonempty square matrix with a dominant eigenvalue separated in modulus.</param>
        /// <param name="iterations">Positive number of iterations; convergence also depends on the starting vector.</param>
        /// <returns>A unit vector V and its Rayleigh quotient D; a zero product returns the current vector and zero.</returns>
        public static (float[] V, float D) Decompose(float[,] matrix, int iterations = 100)
        {
            var d = Iterate(InternalMatrixMath.Copy(matrix, true), iterations);
            return (InternalMatrixMath.Real(d.V), (float)d.D.Real);
        }

        /// <summary>Places an existing power-iteration vector on a diagonal without further iteration.</summary>
        /// <param name="vector">Vector returned by Decompose; these entries are eigenvector components.</param>
        /// <returns>A diagonal matrix containing the vector entries.</returns>
        public static float[,] DiagonalMatrix(float[] vector)
        {
            if (vector == null) throw new ArgumentNullException(nameof(vector));
            return vector.Diag();
        }

        /// <summary>Approximates a dominant right eigenpair by normalized power iteration.</summary>
        /// <param name="matrix">Finite nonempty square matrix with a dominant eigenvalue separated in modulus.</param>
        /// <param name="iterations">Positive number of iterations; convergence also depends on the starting vector.</param>
        /// <returns>A unit vector V and its Rayleigh quotient D; a zero product returns the current vector and zero.</returns>
        public static (Complex32[] V, Complex32 D) Decompose(Complex32[,] matrix, int iterations = 100)
        {
            var d = Iterate(InternalMatrixMath.Copy(matrix, true), iterations);
            return (InternalMatrixMath.Single(d.V), new Complex32((float)d.D.Real, (float)d.D.Imaginary));
        }

        /// <summary>Places an existing power-iteration vector on a diagonal without further iteration.</summary>
        /// <param name="vector">Vector returned by Decompose; these entries are eigenvector components.</param>
        /// <returns>A diagonal matrix containing the vector entries.</returns>
        public static Complex32[,] DiagonalMatrix(Complex32[] vector)
        {
            if (vector == null) throw new ArgumentNullException(nameof(vector));
            return vector.Diag();
        }

        /// <summary>Applies normalized matrix-vector products and computes the final Hermitian Rayleigh quotient.</summary>
        /// <param name="a">Private square matrix.</param>
        /// <param name="iterations">Positive iteration count.</param>
        /// <returns>The unit right vector and corresponding quotient.</returns>
        private static (C[] V, C D) Iterate(C[,] a, int iterations)
        {
            if (iterations < 1) throw new ArgumentOutOfRangeException(nameof(iterations));
            int n = a.GetLength(0);
            var v = new C[n];
            for (int i = 0; i < n; i++) v[i] = 1 / Math.Sqrt(n);
            for (int step = 0; step < iterations; step++)
            {
                var w = InternalMatrixMath.Multiply(a, v);
                double norm = InternalMatrixMath.Norm(w);
                if (norm == 0) return (v, C.Zero);
                InternalMatrixMath.Divide(w, norm);
                v = w;
            }
            C eigenvalue = 0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) eigenvalue += C.Conjugate(v[i]) * a[i, j] * v[j];
            return (v, eigenvalue);
        }
    }
}
