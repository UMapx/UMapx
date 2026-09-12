using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides generalized singular value decomposition for tall matrix pairs</summary>
    public static class GSVD
    {
        /// <summary>Computes A = U1 diag(S1) X and B = U2 diag(S2) X</summary>
        /// <param name="a">Finite m by n matrix with m >= n.</param>
        /// <param name="b">Finite p by n matrix with p >= n. The stacked pair must have full column rank.</param>
        /// <param name="iterations">Positive maximum Jacobi SVD sweeps.</param>
        /// <returns>Orthonormal U1 and U2, nonnegative S1 and S2 with S1^2+S2^2=1, and invertible X.</returns>
        public static (float[,] U1, float[] S1, float[,] U2, float[] S2, float[,] X)
            Decompose(float[,] a, float[,] b, int iterations = 50)
        {
            var d = Factor(InternalMatrixMath.Copy(a), InternalMatrixMath.Copy(b), iterations);
            return (InternalMatrixMath.Real(d.U1), d.S1, InternalMatrixMath.Real(d.U2), d.S2, InternalMatrixMath.Real(d.X));
        }

        /// <summary>Computes A = U1 diag(S1) X and B = U2 diag(S2) X</summary>
        /// <param name="a">Finite m by n matrix with m >= n.</param>
        /// <param name="b">Finite p by n matrix with p >= n. The stacked pair must have full column rank.</param>
        /// <param name="iterations">Positive maximum Jacobi SVD sweeps.</param>
        /// <returns>Orthonormal U1 and U2, nonnegative S1 and S2 with S1^2+S2^2=1, and invertible X.</returns>
        public static (Complex32[,] U1, float[] S1, Complex32[,] U2, float[] S2, Complex32[,] X)
            Decompose(Complex32[,] a, Complex32[,] b, int iterations = 50)
        {
            var d = Factor(InternalMatrixMath.Copy(a), InternalMatrixMath.Copy(b), iterations);
            return (InternalMatrixMath.Single(d.U1), d.S1, InternalMatrixMath.Single(d.U2), d.S2, InternalMatrixMath.Single(d.X));
        }

        /// <summary>Computes generalized singular values from existing diagonal factors</summary>
        /// <param name="s1">Nonnegative first diagonal.</param>
        /// <param name="s2">Nonnegative second diagonal of the same length.</param>
        /// <returns>S1/S2, with positive infinity for positive/zero and NaN for zero/zero.</returns>
        public static float[] GeneralizedSingularValues(float[] s1, float[] s2)
        {
            Validate(s1, s2);
            var result = new float[s1.Length];
            for (int i = 0; i < result.Length; i++) result[i] = s1[i] / s2[i];
            return result;
        }

        /// <summary>Evaluates the GSVD diagonal normalization identity</summary>
        /// <param name="s1">First diagonal.</param>
        /// <param name="s2">Second diagonal of the same length.</param>
        /// <returns>The entries S1[i]^2 + S2[i]^2.</returns>
        public static float[] Identity(float[] s1, float[] s2)
        {
            Validate(s1, s2);
            var result = new float[s1.Length];
            for (int i = 0; i < result.Length; i++) result[i] = (float)((double)s1[i] * s1[i] + (double)s2[i] * s2[i]);
            return result;
        }

        /// <summary>Validates the two real GSVD diagonal factors</summary>
        /// <param name="s1">First nonnegative diagonal.</param>
        /// <param name="s2">Second nonnegative diagonal.</param>
        private static void Validate(float[] s1, float[] s2)
        {
            SVD.Norm(s1); SVD.Norm(s2);
            if (s1.Length != s2.Length) throw new ArgumentException("Diagonal lengths must agree.");
        }

        /// <summary>Combines independently scaled stacked QR with an SVD of the upper orthonormal block</summary>
        /// <param name="a">Private tall first matrix.</param>
        /// <param name="b">Private tall second matrix with the same column count.</param>
        /// <param name="iterations">Jacobi sweep limit.</param>
        /// <returns>Economy GSVD factors, with completed orthonormal columns for zero sine values.</returns>
        private static (C[,] U1, float[] S1, C[,] U2, float[] S2, C[,] X) Factor(C[,] a, C[,] b, int iterations)
        {
            int m = a.GetLength(0), p = b.GetLength(0), n = a.GetLength(1);
            if (b.GetLength(1) != n || m < n || p < n)
                throw new ArgumentException("Both matrices must have the same column count and at least that many rows.");
            // Equalize the input units before QR. Otherwise the smaller block can be lost
            // when the larger block's singular values round to one, leaving its basis unresolved.
            double scaleA = InternalMatrixMath.Max(a), scaleB = InternalMatrixMath.Max(b);
            if (scaleA == 0) scaleA = 1;
            if (scaleB == 0) scaleB = 1;
            var stacked = new C[m + p, n];
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m; i++) stacked[i, j] = a[i, j] / scaleA;
                for (int i = 0; i < p; i++) stacked[m + i, j] = b[i, j] / scaleB;
            }
            var qr = QR.Factor(stacked, full: false);
            var r = InternalMatrixMath.Block(qr.R, n, n);
            double threshold = 32 * InternalMatrixMath.Roundoff * InternalMatrixMath.Max(r);
            for (int i = 0; i < n; i++)
                if (C.Abs(r[i, i]) <= threshold) throw new ArgumentException("The stacked matrix must have full column rank.");
            var svd = SVD.Factor(InternalMatrixMath.Block(qr.Q, m, n), iterations);
            var w = InternalMatrixMath.Multiply(InternalMatrixMath.Block(qr.Q, p, n, m), svd.V);
            var x = InternalMatrixMath.Multiply(InternalMatrixMath.Adjoint(svd.V), r);
            var s1 = new float[n];
            var s2 = new float[n];
            var u2 = new C[p, n];
            var basis = new C[p, n];
            var zero = new bool[n];
            int count = 0;
            for (int j = 0; j < n; j++)
            {
                var column = InternalMatrixMath.Column(w, j);
                double sine = InternalMatrixMath.Norm(column);
                if (sine <= 64 * InternalMatrixMath.Roundoff) { sine = 0; zero[j] = true; }
                else
                {
                    InternalMatrixMath.Divide(column, sine);
                    InternalMatrixMath.SetColumn(u2, j, column);
                    InternalMatrixMath.SetColumn(basis, count, column);
                    count++;
                }
                // Restore each input scale through the diagonal factors and shared X.
                // This preserves S1^2 + S2^2 = 1 without changing either orthonormal basis.
                sine *= scaleB;
                double cosine = svd.S[j] * scaleA;
                double normalization = Math.Sqrt(cosine * cosine + sine * sine);
                s1[j] = (float)(cosine / normalization);
                s2[j] = (float)(sine / normalization);
                for (int k = 0; k < n; k++) x[j, k] *= normalization;
            }
            for (int j = 0; j < n; j++)
                if (zero[j])
                {
                    var column = InternalMatrixMath.Complete(basis, count);
                    InternalMatrixMath.SetColumn(u2, j, column);
                    InternalMatrixMath.SetColumn(basis, count, column);
                    count++;
                }
            return (svd.U, s1, u2, s2, x);
        }
    }
}
