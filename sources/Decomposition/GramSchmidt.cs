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
            var d = Factor(InternalRealMatrixMath.Copy(matrix));
            return (InternalRealMatrixMath.Real(d.Q), InternalRealMatrixMath.Real(d.R));
        }

        /// <summary>Computes an economy QR factorization using reorthogonalized modified Gram-Schmidt.</summary>
        /// <param name="matrix">Finite nonempty matrix with at least as many rows as columns.</param>
        /// <returns>Orthonormal columns Q and upper triangular R, including basis completion for dependent columns.</returns>
        public static (Complex32[,] Q, Complex32[,] R) Decompose(Complex32[,] matrix)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix));
            return (InternalMatrixMath.Single(d.Q), InternalMatrixMath.Single(d.R));
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
                var v = InternalMatrixMath.Column(a, j);
                double original = InternalMatrixMath.Norm(v);
                InternalMatrixMath.Orthogonalize(v, q, j, coefficients: r, column: j);
                double norm = InternalMatrixMath.Norm(v);
                if (norm <= 16 * InternalMatrixMath.Roundoff * original) v = InternalMatrixMath.Complete(q, j);
                else { r[j, j] = norm; InternalMatrixMath.Divide(v, norm); }
                InternalMatrixMath.SetColumn(q, j, v);
            }
            return (q, r);
        }

        /// <summary>Orthogonalizes columns twice and completes the basis at numerical breakdown.</summary>
        /// <param name="a">Private tall or square matrix.</param>
        /// <returns>An economy orthonormal basis and its upper triangular coefficients.</returns>
        private static (double[][] Q, double[][] R) Factor(double[][] a)
        {
            int m = a.Length, n = a[0].Length;
            if (m < n) throw new ArgumentException("Gram-Schmidt requires rows >= columns.");
            var q = InternalRealMatrixMath.Create(m, n);
            var r = InternalRealMatrixMath.Create(n, n);
            for (int j = 0; j < n; j++)
            {
                var v = InternalRealMatrixMath.Column(a, j);
                double original = InternalRealMatrixMath.Norm(v);
                InternalRealMatrixMath.Orthogonalize(v, q, j, coefficients: r, column: j);
                double norm = InternalRealMatrixMath.Norm(v);
                if (norm <= 16 * InternalRealMatrixMath.Roundoff * original) v = InternalRealMatrixMath.Complete(q, j);
                else { r[j][j] = norm; InternalRealMatrixMath.Divide(v, norm); }
                InternalRealMatrixMath.SetColumn(q, j, v);
            }
            return (q, r);
        }
    }
}
