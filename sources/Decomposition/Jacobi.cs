using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines Jacobi decomposition.
    /// </summary>
    /// <remarks>
    /// This class computes the Singular Value Decomposition (SVD) of a real matrix A = U * Σ * Vᵀ
    /// using two algorithmic paths:
    ///   1) If A is (numerically) square and symmetric → symmetric Jacobi eigen-decomposition,
    ///      then converted to SVD (U = Q, Σ = |Λ|, V = Q * diag(sign(Λ))).
    ///   2) Otherwise → one-sided Jacobi SVD (right rotations orthogonalize columns).
    ///
    /// Storage/layout:
    ///   - U is m×n with orthonormal columns (up to numerical tolerance).
    ///   - Σ is length-n (singular values, nonnegative, sorted descending).
    ///   - V is n×n orthonormal.
    ///   
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Singular_value_decomposition
    /// </remarks>
    public class Jacobi
    {
        #region Private data
        private readonly float[][] _U;  // m × n
        private readonly float[] _S;    // n
        private readonly float[][] _V;  // n × n
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes Jacobi decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="eps">Convergence tolerance (0, 1)</param>
        /// <param name="maxSweeps">Maximum number of full sweeps over all index pairs (> 0)</param>
        public Jacobi(float[,] A, float eps = 1e-8f, int maxSweeps = 10)
        {
            // defensive copy
            var X = A.ToJagged();
            int m = X.Length;
            int n = X[0].Length;

            if (m == n && Matrice.IsSymmetric(A))
            {
                // Symmetric path: eigen-Jacobi -> SVD
                var (evals, Q) = JacobiEigenSymmetric(X, eps, maxSweeps); // Q: n×n
                var S = new float[n];
                var U = (float[][])Q.Clone();
                var V = (float[][])Q.Clone();

                for (int j = 0; j < n; j++)
                {
                    float lambda = evals[j];
                    float sigma = Maths.Abs(lambda);
                    S[j] = sigma;
                    float sign = (lambda >= 0f || sigma == 0f) ? 1f : -1f;
                    for (int i = 0; i < n; i++) V[i][j] *= sign;
                }

                SortDescending(ref S, ref U, ref V);
                _U = U; _S = S; _V = V;
            }
            else
            {
                // General path: one-sided Jacobi SVD
                var (U, S, V) = OneSidedJacobiSvd(X, eps, maxSweeps);
                _U = U; _S = S; _V = V;
            }
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets left singular vectors.
        /// </summary>
        public float[,] U => _U.FromJagged();
        /// <summary>
        /// Gets singular values.
        /// </summary>
        public float[] S => _S;
        /// <summary>
        /// Gets right singular vectors.
        /// </summary>
        public float[,] V => _V.FromJagged();
        #endregion

        #region Private voids
        /// <summary>
        /// One-sided Jacobi SVD (float, jagged)
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="eps">Convergence tolerance (0, 1)</param>
        /// <param name="maxSweeps">Maximum number of full sweeps over all index pairs (> 0)</param>
        /// <returns>Result</returns>
        private static (float[][] U, float[] S, float[][] V) OneSidedJacobiSvd(float[][] A, float eps, int maxSweeps)
        {
            int m = A.Length, n = A[0].Length;
            var V = Jagged.Eye(n, n);

            // Frobenius norm squared
            float fro2 = 0f;
            for (int i = 0; i < m; i++)
            {
                var Ai = A[i];
                for (int j = 0; j < n; j++) fro2 += Ai[j] * Ai[j];
            }

            for (int sweep = 0; sweep < maxSweeps; sweep++)
            {
                float corrSum = 0f;

                for (int p = 0; p < n; p++)
                    for (int q = p + 1; q < n; q++)
                    {
                        float alpha = 0f, beta = 0f, gamma = 0f;
                        for (int i = 0; i < m; i++)
                        {
                            float ap = A[i][p], aq = A[i][q];
                            alpha += ap * ap;
                            beta += aq * aq;
                            gamma += ap * aq;
                        }

                        corrSum += gamma * gamma;

                        if (alpha == 0f && beta == 0f) continue;
                        if (Maths.Abs(gamma) <= eps * Maths.Sqrt(alpha * beta)) continue;

                        // robust rotation (float)
                        float tau = (beta - alpha) / (2f * gamma);
                        float at = Maths.Abs(tau);
                        float t = (tau >= 0f ? 1f : -1f) / (at + Maths.Sqrt(1f + at * at));
                        float c = 1f / Maths.Sqrt(1f + t * t);
                        float s = c * t;

                        // A[:,{p,q}] = A * R
                        for (int i = 0; i < m; i++)
                        {
                            float ap = A[i][p], aq = A[i][q];
                            float np = c * ap + s * aq;
                            float nq = -s * ap + c * aq;
                            A[i][p] = np;
                            A[i][q] = nq;
                        }
                        // V = V * R
                        for (int i = 0; i < n; i++)
                        {
                            float vp = V[i][p], vq = V[i][q];
                            float np = c * vp + s * vq;
                            float nq = -s * vp + c * vq;
                            V[i][p] = np;
                            V[i][q] = nq;
                        }
                    }

                if (corrSum <= eps * eps * fro2) break;
            }

            // Build U (m×n) and S (n) from column norms of A
            var Svals = new float[n];
            var U = Jagged.Zero(m, n);
            for (int j = 0; j < n; j++)
            {
                float norm2 = 0f;
                for (int i = 0; i < m; i++) norm2 += A[i][j] * A[i][j];
                float norm = Maths.Sqrt(norm2);
                Svals[j] = norm;
                float inv = (norm > 0f) ? (1f / norm) : 0f;
                for (int i = 0; i < m; i++) U[i][j] = A[i][j] * inv;
            }

            SortDescending(ref Svals, ref U, ref V);
            return (U, Svals, V);
        }
        /// <summary>
        /// Symmetric Jacobi eigen (float, jagged).
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="eps">Convergence tolerance (0, 1)</param>
        /// <param name="maxSweeps">Maximum number of full sweeps over all index pairs (> 0)</param>
        /// <returns>Result</returns>
        private static (float[] evals, float[][] Q) JacobiEigenSymmetric(float[][] A, float eps, int maxSweeps)
        {
            int n = A.Length;
            var Q = Jagged.Eye(n, n);

            for (int sweep = 0; sweep < maxSweeps; sweep++)
            {
                float off2 = 0f;
                for (int i = 0; i < n; i++)
                    for (int j = i + 1; j < n; j++)
                        off2 += A[i][j] * A[i][j];

                if (Maths.Sqrt(off2) <= eps) break;

                for (int i = 0; i < n; i++)
                    for (int j = i + 1; j < n; j++)
                    {
                        float aij = A[i][j];
                        if (Maths.Abs(aij) <= eps) continue;

                        float aii = A[i][i], ajj = A[j][j];

                        float tau = (ajj - aii) / (2f * aij);
                        float at = Maths.Abs(tau);
                        float t = (tau >= 0f ? 1f : -1f) / (at + Maths.Sqrt(1f + at * at));
                        float c = 1f / Maths.Sqrt(1f + t * t);
                        float s = c * t;

                        float aiiNew = aii - t * aij;
                        float ajjNew = ajj + t * aij;
                        A[i][i] = aiiNew;
                        A[j][j] = ajjNew;
                        A[i][j] = 0f;
                        A[j][i] = 0f;

                        for (int k = 0; k < n; k++)
                        {
                            if (k == i || k == j) continue;

                            float aik = A[i][k];
                            float ajk = A[j][k];

                            float aikNew = c * aik - s * ajk;
                            float ajkNew = s * aik + c * ajk;

                            A[i][k] = aikNew;
                            A[k][i] = aikNew; // maintain symmetry in storage
                            A[j][k] = ajkNew;
                            A[k][j] = ajkNew;
                        }

                        // Accumulate Q = Q * G(i,j)
                        for (int k = 0; k < n; k++)
                        {
                            float qki = Q[k][i], qkj = Q[k][j];
                            float qkiN = c * qki - s * qkj;
                            float qkjN = s * qki + c * qkj;
                            Q[k][i] = qkiN;
                            Q[k][j] = qkjN;
                        }
                    }
            }

            var evals = new float[n];
            for (int i = 0; i < n; i++) evals[i] = A[i][i];
            return (evals, Q);
        }
        /// <summary>
        /// Sort Σ desc and permute columns of U (m×n) and V (n×n)
        /// </summary>
        /// <param name="S">Singular values (length n)</param>
        /// <param name="U">Left singular vectors (m×n), columns permuted in place via copy-out</param>
        /// <param name="V">Right singular vectors (n×n), columns permuted accordingly</param>
        private static void SortDescending(ref float[] S, ref float[][] U, ref float[][] V)
        {
            int m = U.Length, n = U[0].Length;
            int[] idx = new int[n];
            for (int i = 0; i < n; i++) idx[i] = i;

            // selection sort by S (desc)
            for (int a = 0; a < n - 1; a++)
            {
                int best = a;
                float bestVal = S[idx[best]];
                for (int b = a + 1; b < n; b++)
                {
                    float val = S[idx[b]];
                    if (val > bestVal) { best = b; bestVal = val; }
                }
                if (best != a) { int t = idx[a]; idx[a] = idx[best]; idx[best] = t; }
            }

            var S2 = new float[n];
            var U2 = Jagged.Zero(m, n);
            var V2 = Jagged.Zero(n, n);

            for (int newj = 0; newj < n; newj++)
            {
                int oldj = idx[newj];
                S2[newj] = S[oldj];

                for (int i = 0; i < m; i++) U2[i][newj] = U[i][oldj];
                for (int i = 0; i < n; i++) V2[i][newj] = V[i][oldj];
            }

            S = S2; U = U2; V = V2;
        }
        #endregion
    }
}
