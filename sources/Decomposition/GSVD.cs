using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Generalized SVD (thin GSVD) for a pair of matrices A (m×n), B (p×n).
    /// Returns orthonormal U (m×n), V (p×n), invertible Q (n×n),
    /// and vectors C[i] = ||A q_i||, S[i] = ||B q_i|| with C[i]^2 + S[i]^2 = 1.
    /// </summary>
    [Serializable]
    public class GSVD
    {
        #region Private data
        private readonly int m, p, n;
        private readonly float eps;

        // thin factors:
        private readonly float[,] u; // m×n, columns orthonormal
        private readonly float[,] v; // p×n, columns orthonormal
        private readonly float[,] q; // n×n, invertible
        private readonly float[] c;  // length n, c_i = ||A q_i||
        private readonly float[] s;  // length n, s_i = ||B q_i||
        #endregion

        #region Public API
        /// <summary>
        /// Left factor with orthonormal columns (m×n).
        /// </summary>
        public float[,] U => u;

        /// <summary>
        /// Right factor with orthonormal columns (p×n).
        /// </summary>
        public float[,] V => v;

        /// <summary>
        /// Common right invertible factor (n×n); columns are C-orthonormal.
        /// </summary>
        public float[,] Q => q;

        /// <summary>
        /// c[i] = ||A q_i||.
        /// </summary>
        public float[] C => c;

        /// <summary>
        /// s[i] = ||B q_i||.
        /// </summary>
        public float[] S => s;
        #endregion

        #region Constructor
        /// <summary>
        /// Computes thin GSVD for A (m×n) and B (p×n).
        /// </summary>
        /// <param name="A">Matrix m×n</param>
        /// <param name="B">Matrix p×n</param>
        /// <param name="eps">Numerical tolerance</param>
        public GSVD(float[,] A, float[,] B, float eps = 1e-8f)
        {
            if (A == null || B == null) throw new ArgumentNullException();
            m = A.GetLength(0);
            n = A.GetLength(1);
            p = B.GetLength(0);
            if (B.GetLength(1) != n) throw new ArgumentException("A and B must have same number of columns");
            if (m <= 0 || n <= 0 || p <= 0) throw new ArgumentException("Invalid dimensions");

            var eps0 = 1e-8f;
            this.eps = eps <= 0 ? eps0 : eps;

            u = new float[m, n];
            v = new float[p, n];
            q = new float[n, n];
            c = new float[n];
            s = new float[n];

            // Core idea:
            // C = A^T A + B^T B (SPD for full column-rank [A;B])
            // Let C = L L^T be Cholesky.
            // S = L^{-T} (A^T A) L^{-1} is symmetric PSD, EVD(S) = Y Λ Y^T.
            // q_i = L^{-1} y_i (C-orthonormalized), then u_i = A q_i / ||A q_i||, v_i = B q_i / ||B q_i||.

            // 1) Build normal matrices:
            var At = Matrice.Transpose(A);
            var Bt = Matrice.Transpose(B);
            var AtA = Matrice.Dot(At, A); // n×n
            var BtB = Matrice.Dot(Bt, B); // n×n
            var Cmat = Matrice.Add(AtA, BtB); // n×n

            // 2) Robust Cholesky with tiny ridge if needed:
            float[,] L = null;
            {
                bool ok = false;
                float ridge = 0f;
                for (int t = 0; t < 5 && !ok; t++)
                {
                    try
                    {
                        var Csym = Symmetrize(t == 0 ? Cmat : AddRidge(Cmat, ridge));
                        var chol = new Cholesky(Csym);
                        L = chol.L; // lower-triangular
                        ok = true;
                    }
                    catch
                    {
                        ridge = (ridge == 0f) ? eps0 : ridge * 10f;
                    }
                }
                if (!ok) throw new ArgumentException("C = A^T A + B^T B is not SPD (even with ridge).");
            }

            // 3) Form S = L^{-T} (A^T A) L^{-1} by solving triangular systems (no explicit inverse):
            var Ssym = new float[n, n];
            {
                var e = new float[n];
                for (int j = 0; j < n; j++)
                {
                    Array.Clear(e, 0, n);
                    e[j] = 1f;

                    var z = new float[n]; SolveLower(L, z, e);          // z = L^{-1} e_j
                    var w = Multiply(AtA, z);                           // w = (A^T A) z
                    var x = new float[n]; SolveUpperT(L, x, w);         // x = L^{-T} w

                    for (int i = 0; i < n; i++) Ssym[i, j] = x[i];
                }
                Ssym = Symmetrize(Ssym);
            }

            // 4) EVD of symmetric S:
            var evd = new EVD(Ssym, this.eps);
            var Y = evd.V;               // n×n, columns y_i (orthonormal)
            var evals = evd.D;           // Complex32[] of length n, ~ real and in [0,1]

            // 5) Sort by λ descending (optional but common):
            var order = ArgsortRealDesc(evals);

            // 6) Build q, u, v, c, s:
            for (int idx = 0; idx < n; idx++)
            {
                int k = order[idx];
                float lambda = evals[k].Real;
                if (lambda < 0f) lambda = 0f;
                if (lambda > 1f) lambda = 1f;

                // q = L^{-1} y, then C-normalize q so that q^T C q = 1
                var y = GetCol(Y, k);
                var qi = new float[n];
                SolveLower(L, qi, y); // L qi = y

                var Cqi = Multiply(Cmat, qi);
                float normC = Maths.Sqrt(Math.Max(Dot(qi, Cqi), 0f));
                if (normC > 0f)
                {
                    float inv = 1f / normC;
                    for (int i = 0; i < n; i++) qi[i] *= inv;
                }
                SetCol(q, idx, qi);

                // columns for U,V and their scales
                var Aq = Multiply(A, qi);
                var Bq = Multiply(B, qi);
                float cVal = Norm2(Aq);
                float sVal = Norm2(Bq);

                // Numerically enforce c^2 + s^2 ≈ 1 (tolerant)
                float sum2 = cVal * cVal + sVal * sVal;
                if (sum2 > 0f)
                {
                    float scale = 1f / Maths.Sqrt(sum2);
                    cVal *= scale;
                    sVal *= scale;
                    for (int i = 0; i < m; i++) Aq[i] *= scale;
                    for (int i = 0; i < p; i++) Bq[i] *= scale;
                }

                c[idx] = cVal;
                s[idx] = sVal;

                float[] uCol = (cVal > this.eps) ? Scale(Aq, 1f / cVal)
                                                 : CreateOrthonormalVector(m, u, idx, this.eps);
                float[] vCol = (sVal > this.eps) ? Scale(Bq, 1f / sVal)
                                                 : CreateOrthonormalVector(p, v, idx, this.eps);

                // Light modified Gram–Schmidt re-orthogonalization (helps when s ~ 0)
                if (cVal <= this.eps) MgsReorth(this.u, idx, uCol, eps0);
                if (sVal <= this.eps) MgsReorth(this.v, idx, vCol, eps0);

                SetCol(u, idx, uCol);
                SetCol(v, idx, vCol);
            }
        }
        #endregion

        #region Helpers (linear algebra)
        /// <summary>
        /// Symmetrizes a square matrix: 0.5*(M + M^T).
        /// </summary>
        private static float[,] Symmetrize(float[,] M)
        {
            int n = M.GetLength(0), m = M.GetLength(1);
            if (n != m) throw new ArgumentException("Symmetrize expects square matrix");
            var S = new float[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    float v = 0.5f * (M[i, j] + M[j, i]);
                    S[i, j] = v;
                    S[j, i] = v;
                }
            }
            return S;
        }

        /// <summary>
        /// Adds ridge (tau * I) to a square matrix.
        /// </summary>
        private static float[,] AddRidge(float[,] M, float tau)
        {
            int n = M.GetLength(0), m = M.GetLength(1);
            if (n != m) throw new ArgumentException("AddRidge expects square matrix");
            var R = (float[,])M.Clone();
            for (int i = 0; i < n; i++) R[i, i] += tau;
            return R;
        }

        /// <summary>
        /// Forward substitution: solve L x = b for lower-triangular L (n×n).
        /// </summary>
        private static void SolveLower(float[,] L, float[] x, float[] b)
        {
            int n = L.GetLength(0);
            for (int i = 0; i < n; i++)
            {
                double sum = b[i];
                for (int k = 0; k < i; k++) sum -= L[i, k] * x[k];
                x[i] = (float)(sum / L[i, i]);
            }
        }

        /// <summary>
        /// Back substitution for transpose: solve L^T x = b with L lower-triangular.
        /// </summary>
        private static void SolveUpperT(float[,] L, float[] x, float[] b)
        {
            int n = L.GetLength(0);
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = b[i];
                for (int k = i + 1; k < n; k++) sum -= L[k, i] * x[k];
                x[i] = (float)(sum / L[i, i]);
            }
        }

        /// <summary>
        /// Matrix-vector product (M v).
        /// </summary>
        private static float[] Multiply(float[,] M, float[] v)
        {
            int r = M.GetLength(0), c = M.GetLength(1);
            if (v.Length != c) throw new ArgumentException();
            var y = new float[r];

            for (int i = 0; i < r; i++)
            {
                double s = 0.0;
                for (int j = 0; j < c; j++) s += M[i, j] * v[j];
                y[i] = (float)s;
            }
            return y;
        }

        /// <summary>
        /// Euclidean 2-norm of vector.
        /// </summary>
        private static float Norm2(float[] v)
        {
            double s = 0.0;
            for (int i = 0; i < v.Length; i++) s += (double)v[i] * v[i];
            return (float)Math.Sqrt(s);
        }

        /// <summary>
        /// Dot product of vectors.
        /// </summary>
        private static float Dot(float[] a, float[] b)
        {
            if (a.Length != b.Length) throw new ArgumentException();
            double s = 0.0;
            for (int i = 0; i < a.Length; i++) s += (double)a[i] * b[i];
            return (float)s;
        }

        /// <summary>
        /// Returns a scaled copy of vector.
        /// </summary>
        private static float[] Scale(float[] v, float alpha)
        {
            var w = new float[v.Length];
            for (int i = 0; i < v.Length; i++) w[i] = alpha * v[i];
            return w;
        }

        /// <summary>
        /// Gets a column j from matrix M as a vector.
        /// </summary>
        private static float[] GetCol(float[,] M, int j)
        {
            int r = M.GetLength(0);
            var v = new float[r];
            for (int i = 0; i < r; i++) v[i] = M[i, j];
            return v;
        }

        /// <summary>
        /// Sets column j of matrix M to vector v.
        /// </summary>
        private static void SetCol(float[,] M, int j, float[] v)
        {
            int r = M.GetLength(0);
            if (v.Length != r) throw new ArgumentException();
            for (int i = 0; i < r; i++) M[i, j] = v[i];
        }

        /// <summary>
        /// Modified Gram–Schmidt re-orthogonalization of vector 'v' against first 'uptoCol' columns of M.
        /// </summary>
        private static void MgsReorth(float[,] M, int uptoCol, float[] v, float tol)
        {
            for (int j = 0; j < uptoCol; j++)
            {
                var w = GetCol(M, j);
                float alpha = Dot(w, v);
                for (int i = 0; i < v.Length; i++) v[i] -= alpha * w[i];
            }
            float nrm = Norm2(v);
            if (nrm > tol)
            {
                float inv = 1f / nrm;
                for (int i = 0; i < v.Length; i++) v[i] *= inv;
            }
        }

        /// <summary>
        /// Creates an orthonormal vector completing the current set of columns in M (fallback when scale ~ 0).
        /// </summary>
        private static float[] CreateOrthonormalVector(int dim, float[,] M, int uptoCol, float tol)
        {
            var rand = new Random(1234 + uptoCol);
            var v = new float[dim];
            for (int i = 0; i < dim; i++) v[i] = (float)(rand.NextDouble() - 0.5);
            MgsReorth(M, uptoCol, v, tol);
            // If degenerate, pick a canonical basis direction not spanned yet:
            if (Norm2(v) < tol)
            {
                for (int k = 0; k < dim; k++)
                {
                    Array.Clear(v, 0, dim);
                    v[k] = 1f;
                    MgsReorth(M, uptoCol, v, tol);
                    if (Norm2(v) >= tol) break;
                }
            }
            return v;
        }

        /// <summary>
        /// Returns indices that sort eigenvalues (Complex32) by descending real part.
        /// </summary>
        private static int[] ArgsortRealDesc(Complex32[] vals)
        {
            int n = vals.Length;
            var idx = new int[n];
            for (int i = 0; i < n; i++) idx[i] = i;
            Array.Sort(idx, (i, j) => -vals[i].Real.CompareTo(vals[j].Real));
            return idx;
        }
        #endregion
    }
}
