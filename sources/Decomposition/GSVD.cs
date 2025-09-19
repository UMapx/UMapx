using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines the generalized singular value decomposition (GSVD).
    /// </summary>
    /// <remarks>
    /// Based on the algorithm proposed by Paige and Saunders and implemented in
    /// https://github.com/baptistefraisse/gsvd. The decomposition factorizes two
    /// real matrices A and B with the same column dimension according to
    /// A = U1 * S1 * X and B = U2 * S2 * X, where U1 and U2 have orthonormal
    /// columns, S1 and S2 are diagonal (with non-negative entries) and X is an
    /// invertible matrix.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Generalized_singular_value_decomposition
    /// </remarks>
    [Serializable]
    public class GSVD
    {
        #region Private data
        private readonly float[,] u1;
        private readonly float[,] u2;
        private readonly float[] s1;
        private readonly float[] s2;
        private readonly float[,] x;
        private readonly float[] gamma;
        private readonly float tolerance;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes generalized singular value decomposition.
        /// </summary>
        /// <param name="A">Matrix A (m x n, with m ≥ n)</param>
        /// <param name="B">Matrix B (p x n, with p ≥ n)</param>
        /// <param name="iterations">Number of iterations of SVD</param>
        public GSVD(float[,] A, float[,] B, int iterations = 10)
        {
            int aRows = A.GetLength(0);
            int aCols = A.GetLength(1);
            int bRows = B.GetLength(0);
            int bCols = B.GetLength(1);

            if (aRows < 2 || aCols < 2 || bRows < 2 || bCols < 2)
                throw new ArgumentException("Error: A and B must size > 1");

            if (aCols != bCols)
                throw new ArgumentException("Error: A and B must have matching column ranks");

            if (aRows < aCols || bRows < bCols)
                throw new ArgumentException("Error: A and B must have full column rank (m>n)");

            // Step 1: QR decomposition of [A; B]
            float[,] stacked = A.Concat(B, Direction.Vertical);
            var qr = new QR(stacked);
            float[,] Q = qr.Q; // (aRows + bRows) x columns
            float[,] R = qr.R; // columns x columns (upper triangular)

            float[,] q1 = ExtractRows(Q, 0, aRows);
            float[,] q2 = ExtractRows(Q, aRows, bRows);

            // Step 2: SVD of the upper part of Q
            var svd = new SVD(q1, iterations); 
            int columns = aCols;
            float[,] U1 = svd.U;
            float[,] V = svd.V;
            float[] singular = svd.S;
            float[] S1 = singular;
            float[] sigma2 = Matrice.Zero(columns);

            for (int i = 0; i < columns; i++)
            {
                float value = 1.0f - singular[i] * singular[i];
                sigma2[i] = value > 0 ? Maths.Sqrt(value) : 0.0f;
            }

            // Step 3: compute U2 = Q2 * V * inv(S2)
            float[] S2 = sigma2;
            float[,] U2;
            float[,] temp = V.Dot(S2, true);
            U2 = q2.Dot(temp);

            // Step 4: diagonal matrix H built from Vᵀ * R
            float[,] VT = V.Transpose();
            float[,] VTR = VT.Dot(R);
            float[,] Hproduct = VTR.Dot(VTR.Transpose());
            float[] hDiag = Hproduct.Diag();
            float[] sqrtH = Matrice.Compute(hDiag, Maths.Sqrt);
            float[,] sqrtHDiag = sqrtH.Diag();
            float[,] invSqrtHDiag = sqrtH.Invert().Diag();

            // Step 5: refine S1, S2 and compute X
            S1 = S1.Dot(sqrtHDiag);
            S2 = S2.Dot(sqrtHDiag);
            float[,] X = invSqrtHDiag.Dot(VTR);

            // Step 6: refine X to satisfy S1^2 + S2^2 = 1
            int n = S1.Length;
            for (int i = 0; i < n; i++)
            {
                float d = Maths.Hypotenuse(S1[i], S2[i]);

                if (d > 0f)
                {
                    float inv = 1f / d;
                    S1[i] *= inv;
                    S2[i] *= inv;

                    // X' = D * X  ⇒  scale i-th row of X by d
                    for (int j = 0; j < n; j++)
                        X[i, j] *= d;
                }
            }

            // Step 7: generalized singular values gamma = diag(S1)/diag(S2)
            float[] Gamma = Matrice.Zero(columns);

            for (int i = 0; i < columns; i++)
            {
                float denom = S2[i];
                Gamma[i] = Maths.Abs(denom) <= 0 ? 0.0f : S1[i] / denom;
            }

            // Store results
            this.u1 = U1;
            this.u2 = U2;
            this.s1 = S1;
            this.s2 = S2;
            this.x = X;
            this.gamma = Gamma;
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Returns matrix U1.
        /// </summary>
        public float[,] U1 { get { return u1; } }
        /// <summary>
        /// Returns matrix U2.
        /// </summary>
        public float[,] U2 { get { return u2; } }
        /// <summary>
        /// Returns diagonal matrix S1.
        /// </summary>
        public float[] S1 { get { return s1; } }
        /// <summary>
        /// Returns diagonal matrix S2.
        /// </summary>
        public float[] S2 { get { return s2; } }
        /// <summary>
        /// Returns matrix X.
        /// </summary>
        public float[,] X { get { return x; } }
        /// <summary>
        /// Returns the generalized singular values.
        /// </summary>
        public float[] Gamma { get { return gamma; } }
        /// <summary>
        /// Returns the identity vector.
        /// </summary>
        public float[] Identity
        {
            get
            {
                var length = s1.Length;
                var one = new float[length];
                for (int i = 0; i < length; i++) one[i] = s1[i] * s1[i] + s2[i] * s2[i];
                return one;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Extracts a consecutive block of rows from a matrix.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="rowStart">Row start</param>
        /// <param name="rowCount">Row count</param>
        /// <returns>Matrix</returns>
        private static float[,] ExtractRows(float[,] matrix, int rowStart, int rowCount)
        {
            int cols = matrix.GetLength(1);
            float[,] result = new float[rowCount, cols];

            for (int i = 0; i < rowCount; i++)
            {
                int row = rowStart + i;
                for (int j = 0; j < cols; j++)
                {
                    result[i, j] = matrix[row, j];
                }
            }

            return result;
        }
        #endregion
    }
}
