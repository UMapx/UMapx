using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Fast Chebyshev transform of the first kind (Type-I), orthonormal.
    /// </summary>
    /// <remarks>
    /// Uses one FFT of the even-symmetric extension (length 2M, where M=N-1) and exact endpoint handling
    /// to reproduce the same result as <see cref="ChebyshevTransform"/>. Forward and inverse are numerically
    /// stable and bitwise consistent with the matrix version up to floating-point rounding.
    /// </remarks>
    [Serializable]
    public class FastChebyshevTransform : TransformBaseFloat, ITransform
    {
        #region Private data
        private readonly FastFourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Chebyshev transform (Type-I, orthonormal).
        /// </summary>
        /// <param name="direction">Processing direction (unused for 1D)</param>
        public FastChebyshevTransform(Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(false, Direction.Both);
            this.Direction = direction;
        }
        #endregion

        #region Fast Chebyshev Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override float[] Forward(float[] A)
        {
            int N = A.Length;
            if (N == 0) return Array.Empty<float>();
            if (N == 1) return new float[] { A[0] };

            int M = N - 1;               // DCT-I size parameter
            int L = 2 * M;               // FFT length (even extension)
            float invSqrt2 = 1.0f / Maths.Sqrt2;

            // Build even-symmetric extension z[0..2M-1] = [x0, x1, ..., xM, xM-1, ..., x1]
            var z = new Complex32[L];
            z[0] = A[0];
            for (int j = 1; j < M; j++) z[j] = A[j];
            z[M] = A[M];
            for (int j = M + 1; j < L; j++) z[j] = A[L - j];

            var Y = FFT.Forward(z); // Unscaled DFT

            float[] c = new float[N];
            float rowScale, Ek, Zk, t;
            float root = Maths.Sqrt(2.0f / M); // √(2/M)

            for (int k = 0; k <= M; k++)
            {
                Zk = Y[k].Real;                                      // x0 + (-1)^k xM + 2∑_{j=1}^{M-1} xj cos(πjk/M)
                Ek = A[0] + (((k & 1) == 0) ? A[M] : -A[M]);         // x0 + (-1)^k xM
                t = 0.5f * Zk + (invSqrt2 - 0.5f) * Ek;              // ∑ a(j)*xj*cos(πjk/M)
                rowScale = ((k == 0) || (k == M)) ? invSqrt2 : 1.0f; // a(k)
                c[k] = root * rowScale * t;                          // c_k = √(2/M)*a(k)*∑ a(j)xj cos
            }
            return c;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override float[] Backward(float[] B)
        {
            // Orthonormal Chebyshev (DCT-I) matrix is symmetric and orthogonal:
            // H^T = H and H^T H = I => H^{-1} = H. Therefore the inverse equals the forward.
            // This guarantees exact equality with the matrix-based Backward (up to FP rounding).
            return Forward(B);
        }
        #endregion
    }
}
