using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines generalized singular value decomposition for a pair of real matrices.
    /// </summary>
    /// <remarks>
    /// Finds orthogonal matrices U and V, and a nonsingular matrix Q such that
    /// Uᵀ * A * Q = diag(C) and Vᵀ * B * Q = diag(S) with non-negative diagonals
    /// satisfying Cᵀ * C + Sᵀ * S = I. Only single precision (float) arithmetic is supported.
    /// </remarks>
    [Serializable]
    public class GSVD
    {
        #region Private data
        private readonly int m;
        private readonly int p;
        private readonly int n;
        private readonly float[,] u;
        private readonly float[,] v;
        private readonly float[,] q;
        private readonly float[] c;
        private readonly float[] s;
        private readonly float eps;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes generalized singular value decomposition.
        /// </summary>
        /// <param name="A">Matrix A (m × n)</param>
        /// <param name="B">Matrix B (p × n)</param>
        /// <param name="eps">Tolerance used when detecting rank deficiencies</param>
        public GSVD(float[,] A, float[,] B, float eps = 1e-8f)
        {
            if (A == null || B == null)
                throw new ArgumentNullException(A == null ? nameof(A) : nameof(B));

            if (A.GetLength(1) != B.GetLength(1))
                throw new ArgumentException("Matrices must have the same number of columns");

            this.m = A.GetLength(0);
            this.p = B.GetLength(0);
            this.n = A.GetLength(1);
            this.eps = (eps < 0) ? 0f : eps;

            // Form the symmetric positive definite matrix C = AᵀA + BᵀB.
            float[,] ata = Matrice.Dot(A.Transpose(), A);
            float[,] btb = Matrice.Dot(B.Transpose(), B);
            float[,] cMat = ata.Add(btb);

            // Cholesky factorization of C.
            var chol = new Cholesky(cMat);
            float[,] l = chol.L;
            float[,] lInv = l.Invert();
            float[,] lInvT = lInv.Transpose();

            // Reduce the generalized eigenvalue problem to a standard symmetric eigenproblem.
            float[,] tmp = Matrice.Dot(ata, lInv);
            float[,] sym = Symmetrize(Matrice.Dot(lInvT, tmp));

            // Solve the symmetric eigenproblem.
            var evd = new EVD(sym, this.eps);
            float[,] eigenvectors = evd.V;
            Complex32[] eigenvalues = evd.D;

            this.u = new float[m, n];
            this.v = new float[p, n];
            this.q = new float[n, n];
            this.c = new float[n];
            this.s = new float[n];

            int[] order = new int[n];
            for (int i = 0; i < n; i++)
            {
                order[i] = i;
            }
            float[] eigenReal = new float[n];
            for (int i = 0; i < n; i++)
            {
                eigenReal[i] = eigenvalues[i].Real;
            }
            Array.Sort(order, (a, b) => eigenReal[b].CompareTo(eigenReal[a]));

            for (int idx = 0; idx < n; idx++)
            {
                int k = order[idx];
                Complex32 lambda = eigenvalues[k];
                if (Math.Abs(lambda.Imag) > 1e-3f)
                    throw new InvalidOperationException("Generalized singular values are complex");

                float[] y = eigenvectors.GetCol(k);
                float[] qVec = Multiply(lInv, y);

                float[] cTimesQ = Multiply(cMat, qVec);
                float normC = Maths.Sqrt(Math.Max(qVec.Dot(cTimesQ), 0f));
                if (normC < this.eps)
                    throw new InvalidOperationException("Failed to normalize GSVD basis vector");

                float invNormC = 1f / normC;
                float[] qNormalized = Scale(qVec, invNormC);
                SetColumn(this.q, idx, qNormalized);

                float[] aCol = Multiply(A, qNormalized);
                float[] bCol = Multiply(B, qNormalized);

                float cVal = Matrice.Norm(aCol);
                float sVal = Matrice.Norm(bCol);

                if (cVal <= this.eps) cVal = 0f;
                if (sVal <= this.eps) sVal = 0f;

                this.c[idx] = cVal;
                this.s[idx] = sVal;

                float[] uCol = (cVal > this.eps) ? Scale(aCol, 1f / cVal) : CreateOrthonormalVector(m, this.u, idx, this.eps);
                float[] vCol = (sVal > this.eps) ? Scale(bCol, 1f / sVal) : CreateOrthonormalVector(p, this.v, idx, this.eps);

                SetColumn(this.u, idx, uCol);
                SetColumn(this.v, idx, vCol);
            }
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Gets the left orthogonal factor U (m × n with orthonormal columns).
        /// </summary>
        public float[,] U => u;
        /// <summary>
        /// Gets the right orthogonal factor V (p × n with orthonormal columns).
        /// </summary>
        public float[,] V => v;
        /// <summary>
        /// Gets the nonsingular transformation matrix Q.
        /// </summary>
        public float[,] Q => q;
        /// <summary>
        /// Gets the generalized cosines (diagonal of C).
        /// </summary>
        public float[] C => c;
        /// <summary>
        /// Gets the generalized sines (diagonal of S).
        /// </summary>
        public float[] S => s;
        #endregion

        #region Private voids
        private static float[,] Symmetrize(float[,] matrix)
        {
            int n = matrix.GetLength(0);
            float[,] result = new float[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    float value = 0.5f * (matrix[i, j] + matrix[j, i]);
                    result[i, j] = value;
                    result[j, i] = value;
                }
            }
            return result;
        }

        private static float[] Multiply(float[,] matrix, float[] vector)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            float[] result = new float[rows];
            for (int i = 0; i < rows; i++)
            {
                float sum = 0f;
                for (int j = 0; j < cols; j++)
                {
                    sum += matrix[i, j] * vector[j];
                }
                result[i] = sum;
            }
            return result;
        }

        private static float[] Scale(float[] vector, float scale)
        {
            float[] result = new float[vector.Length];
            for (int i = 0; i < vector.Length; i++)
            {
                result[i] = vector[i] * scale;
            }
            return result;
        }

        private static void SetColumn(float[,] matrix, int column, float[] data)
        {
            int rows = matrix.GetLength(0);
            for (int i = 0; i < rows; i++)
            {
                matrix[i, column] = (i < data.Length) ? data[i] : 0f;
            }
        }

        private static float[] CreateOrthonormalVector(int size, float[,] existing, int count, float eps)
        {
            float[] candidate = new float[size];
            for (int attempt = 0; attempt < size; attempt++)
            {
                for (int i = 0; i < size; i++)
                    candidate[i] = 0f;
                candidate[attempt] = 1f;

                for (int j = 0; j < count; j++)
                {
                    float[] prev = existing.GetCol(j);
                    float proj = prev.Dot(candidate);
                    for (int i = 0; i < size; i++)
                        candidate[i] -= proj * prev[i];
                }

                float norm = candidate.Norm();
                if (norm > eps)
                {
                    for (int i = 0; i < size; i++)
                        candidate[i] /= norm;
                    return candidate;
                }
            }

            for (int i = 0; i < size; i++)
                candidate[i] = 0f;
            candidate[0] = 1f;

            for (int j = 0; j < count; j++)
            {
                float[] prev = existing.GetCol(j);
                float proj = prev.Dot(candidate);
                for (int i = 0; i < size; i++)
                    candidate[i] -= proj * prev[i];
            }

            float fallbackNorm = candidate.Norm();
            if (fallbackNorm > 0f)
            {
                for (int i = 0; i < size; i++)
                    candidate[i] /= fallbackNorm;
            }

            return candidate;
        }

        #endregion
    }
}
