using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides real and complex Arnoldi reduction.</summary>
    public static class Arnoldi
    {
        /// <summary>Computes a full Arnoldi reduction A = Q H Q^T.</summary>
        /// <param name="matrix">Finite nonempty square input, not modified.</param>
        /// <returns>An orthonormal basis Q and upper Hessenberg H; invariant-subspace breakdown starts a new orthogonal block.</returns>
        public static (float[,] Q, float[,] H) Decompose(float[,] matrix)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix, true), true);
            return (InternalMatrixMath.Real(d.Q), InternalMatrixMath.Real(d.H));
        }

        /// <summary>Computes a full Arnoldi reduction A = Q H Q^H.</summary>
        /// <param name="matrix">Finite nonempty square input, not modified.</param>
        /// <returns>An orthonormal basis Q and upper Hessenberg H; invariant-subspace breakdown starts a new orthogonal block.</returns>
        public static (Complex32[,] Q, Complex32[,] H) Decompose(Complex32[,] matrix)
        {
            var d = Factor(InternalMatrixMath.Copy(matrix, true), true);
            return (InternalMatrixMath.Single(d.Q), InternalMatrixMath.Single(d.H));
        }

        /// <summary>Builds a complete Krylov basis with deterministic orthogonal restarts.</summary>
        /// <param name="a">Private square work matrix.</param>
        /// <param name="full">Whether to use a second reorthogonalization pass.</param>
        /// <returns>Q and H satisfying A Q = Q H to working precision.</returns>
        internal static (C[,] Q, C[,] H) Factor(C[,] a, bool full)
        {
            int n = a.GetLength(0);
            var q = new C[n, n];
            var h = new C[n, n];
            for (int i = 0; i < n; i++) q[i, 0] = 1 / Math.Sqrt(n);
            double threshold = 32 * InternalMatrixMath.Roundoff * n * InternalMatrixMath.Max(a);
            for (int k = 0; k < n; k++)
            {
                var v = new C[n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) v[i] += a[i, j] * q[j, k];
                InternalMatrixMath.Orthogonalize(v, q, k + 1, full ? 2 : 1, h, k);
                if (k + 1 == n) continue;
                double norm = InternalMatrixMath.Norm(v);
                if (norm <= threshold) v = InternalMatrixMath.Complete(q, k + 1);
                else { h[k + 1, k] = norm; InternalMatrixMath.Divide(v, norm); }
                InternalMatrixMath.SetColumn(q, k + 1, v);
            }
            return (q, h);
        }
    }
}
