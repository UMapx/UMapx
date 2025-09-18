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
        /// <param name="eps">Tolerance used to treat near-zero values</param>
        public GSVD(float[,] A, float[,] B, float eps = 1e-8f)
        {
            if (A == null)
                throw new ArgumentNullException(nameof(A));
            if (B == null)
                throw new ArgumentNullException(nameof(B));
            if (eps <= 0)
                throw new ArgumentOutOfRangeException(nameof(eps));

            tolerance = eps;

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
            SVD svd = new SVD(q1);
            float[,] U1 = svd.U;
            float[,] V = svd.V;
            float[] singular = svd.S;

            int columns = aCols;
            float[] S1 = singular;
            float[] sigma2 = new float[columns];

            for (int i = 0; i < columns; i++)
            {
                float value = 1.0f - singular[i] * singular[i];
                if (value < -tolerance)
                    throw new ArgumentException("Error: A and B similarly singular");

                sigma2[i] = value > tolerance ? Maths.Sqrt(value) : 0.0f;
            }

            // Step 3: compute U2 = Q2 * V * inv(S2)
            float[] S2 = sigma2;
            float[] S2inv = InvertDiagonal(sigma2);
            float[,] U2;
            try
            {
                float[,] temp = V.Dot(S2inv);
                U2 = q2.Dot(temp);
            }
            catch (Exception)
            {
                throw new ArgumentException("Error: A and B similarly singular");
            }

            // Step 4: diagonal matrix H built from Vᵀ * R
            float[,] VT = V.Transpose();
            float[,] VTR = VT.Dot(R);
            float[,] Hproduct = VTR.Dot(VTR.Transpose());
            float[] hDiag = Hproduct.Diag();

            float[,] sqrtH = Matrice.Zero(columns, columns);
            float[,] invSqrtH = Matrice.Zero(columns, columns);

            for (int i = 0; i < columns; i++)
            {
                float value = hDiag[i];
                if (value < 0 && Maths.Abs(value) <= tolerance)
                {
                    value = 0.0f;
                }

                if (value < -tolerance)
                    throw new ArgumentException("Error: A and B similarly singular");

                float sqrtVal = value > tolerance ? Maths.Sqrt(value) : 0.0f;
                sqrtH[i, i] = sqrtVal;
                invSqrtH[i, i] = sqrtVal > tolerance ? 1.0f / sqrtVal : 0.0f;
            }

            // Step 5: refine S1, S2 and compute X
            S1 = S1.Dot(sqrtH);
            S2 = S2.Dot(sqrtH);
            float[,] X = invSqrtH.Dot(VTR);

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
            float[] Gamma = new float[columns];

            for (int i = 0; i < columns; i++)
            {
                float denom = S2[i];
                Gamma[i] = Maths.Abs(denom) <= tolerance ? 0.0f : S1[i] / denom;
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
        /// <summary>
        /// Returns the inverse of a diagonal matrix represented by its diagonal.
        /// </summary>
        /// <param name="diagonal">Array</param>
        /// <returns>Array</returns>
        private float[] InvertDiagonal(float[] diagonal)
        {
            int n = diagonal.Length;
            float[] result = Matrice.Zero(n);

            for (int i = 0; i < n; i++)
            {
                float value = diagonal[i];
                if (Maths.Abs(value) <= tolerance)
                    throw new ArgumentException("Error: A and B similarly singular");

                result[i] = 1.0f / value;
            }

            return result;
        }
        #endregion
    }
}
