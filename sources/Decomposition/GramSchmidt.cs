using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides modified Gram-Schmidt orthogonalization.</summary>
    public static class GramSchmidt
    {
        /// <summary>Computes an economy QR factorization using reorthogonalized modified Gram-Schmidt.</summary>
        /// <param name="matrix">Finite nonempty matrix with at least as many rows as columns.</param>
        /// <returns>Orthonormal columns Q and upper triangular R, including basis completion for dependent columns.</returns>
        public static (float[,] Q, float[,] R) Decompose(float[,] matrix)
        {
            var d = Factor(MatrixMath.Copy(matrix));
            return (MatrixMath.Real(d.Q), MatrixMath.Real(d.R));
        }

        /// <summary>Computes an economy QR factorization using reorthogonalized modified Gram-Schmidt.</summary>
        /// <param name="matrix">Finite nonempty matrix with at least as many rows as columns.</param>
        /// <returns>Orthonormal columns Q and upper triangular R, including basis completion for dependent columns.</returns>
        public static (Complex32[,] Q, Complex32[,] R) Decompose(Complex32[,] matrix)
        {
            var d = Factor(MatrixMath.Copy(matrix));
            return (MatrixMath.Single(d.Q), MatrixMath.Single(d.R));
        }

        /// <summary>Orthogonalizes columns twice and completes the basis at numerical breakdown.</summary>
        /// <param name="a">Private tall or square matrix.</param>
        /// <returns>An economy orthonormal basis and its upper triangular coefficients.</returns>
        private static (C[,] Q, C[,] R) Factor(C[,] a)
        {
            int m = a.GetLength(0), n = a.GetLength(1);
            if (m < n) throw new ArgumentException("Gram-Schmidt requires rows >= columns.");
            var q = new C[m, n];
            var r = new C[n, n];
            for (int j = 0; j < n; j++)
            {
                var v = new C[m];
                for (int i = 0; i < m; i++) v[i] = a[i, j];
                double original = MatrixMath.Norm(v);
                for (int pass = 0; pass < 2; pass++)
                    for (int k = 0; k < j; k++)
                    {
                        C dot = 0;
                        for (int i = 0; i < m; i++) dot += C.Conjugate(q[i, k]) * v[i];
                        r[k, j] += dot;
                        for (int i = 0; i < m; i++) v[i] -= q[i, k] * dot;
                    }
                double norm = MatrixMath.Norm(v);
                if (norm <= 16 * MatrixMath.Roundoff * original) v = MatrixMath.Complete(q, j);
                else { r[j, j] = norm; for (int i = 0; i < m; i++) v[i] /= norm; }
                for (int i = 0; i < m; i++) q[i, j] = v[i];
            }
            return (q, r);
        }
    }
}
