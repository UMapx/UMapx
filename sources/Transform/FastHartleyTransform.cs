using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Pure fast Hartley transform (FHT), no FFT coupling.
    /// Radix-2 DIT recursion with exact Hartley twiddle algebra.
    /// </summary>
    /// <remarks>
    /// Core combine (for N = 2M):
    ///   X[k]     = E[k] + cos(πk/M) * O[k] + sin(πk/M) * O[(M - k) % M]
    ///   X[k + M] = E[k] - cos(πk/M) * O[k] - sin(πk/M) * O[(M - k) % M]
    /// This follows the classic Bracewell/Sorensen FHT identities (self-inverse up to scaling).
    /// </remarks>
    [Serializable]
    public class FastHartleyTransform : TransformBaseFloat, ITransform
    {
        #region Private data
        private bool normalized;
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Hartley transform (pure FHT, no FFT).
        /// </summary>
        /// <param name="normalized">Apply 1/sqrt(N) scaling (orthonormal form)</param>
        /// <param name="direction">Processing direction (kept for API parity)</param>
        public FastHartleyTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.normalized = normalized;
            this.direction = direction;
        }

        /// <summary>
        /// Normalized transform or not (1/sqrt(N) scaling on both forward/backward).
        /// </summary>
        public bool Normalized
        {
            get => normalized;
            set => normalized = value;
        }

        /// <summary>
        /// Gets or sets the processing direction (placeholder for matrix paths).
        /// </summary>
        public override Direction Direction
        {
            get => direction;
            set => direction = value;
        }
        #endregion

        #region Fast Hartley Transform (1D, float)
        /// <summary>
        /// Forward Hartley transform (real → real), O(N log N) for even N via pure FHT.
        /// For odd N (or tiny N), falls back to safe O(N^2) Hartley kernel.
        /// </summary>
        public override float[] Forward(float[] A)
        {
            if (A == null) return null;
            int N = A.Length;
            if (N == 0) return Array.Empty<float>();

            var outv = new float[N];
            FHT(A, 0, 1, N, outv, 0);

            if (normalized)
            {
                float s = 1f / Maths.Sqrt(N);
                for (int i = 0; i < N; i++) outv[i] *= s;
            }
            return outv;
        }

        /// <summary>
        /// Backward Hartley transform (real → real), identical to forward up to scaling.
        /// </summary>
        public override float[] Backward(float[] B)
        {
            if (B == null) return null;
            int N = B.Length;
            if (N == 0) return Array.Empty<float>();

            var outv = new float[N];
            FHT(B, 0, 1, N, outv, 0);

            if (normalized)
            {
                float s = 1f / Maths.Sqrt(N);
                for (int i = 0; i < N; i++) outv[i] *= s;
            }
            return outv;
        }
        #endregion

        #region Core FHT
        /// <summary>
        /// Radix-2 DIT FHT for even n; safe O(n^2) direct Hartley for odd n.
        /// Writes n outputs contiguously to dst[dOff .. dOff + n - 1].
        /// </summary>
        private static void FHT(float[] src, int sOff, int sStride, int n, float[] dst, int dOff)
        {
            if (n <= 0) return;

            if (n == 1)
            {
                dst[dOff] = src[sOff];
                return;
            }

            if ((n & 1) == 0)
            {
                // n = 2M: compute FHT on evens and odds (length M each).
                int M = n >> 1;
                int dMid = dOff + M;

                // E ← FHT(evens), O ← FHT(odds)
                FHT(src, sOff, sStride << 1, M, dst, dOff);
                FHT(src, sOff + sStride, sStride << 1, M, dst, dMid);

                // k = 0 (degenerate pair)
                {
                    float E0 = dst[dOff + 0];
                    float O0 = dst[dMid + 0];
                    dst[dOff + 0] = E0 + O0; // X[0]
                    dst[dOff + M + 0] = E0 - O0; // X[M]
                }

                int half = M >> 1;

                // 1 .. half-1 : pair k with (M - k)
                for (int k = 1; k < half; k++)
                {
                    int km = M - k;

                    // Twiddles β = π k / M
                    double beta = Math.PI * k / M;
                    float c = (float)Math.Cos(beta);
                    float s = (float)Math.Sin(beta);

                    // Read inputs first (avoid overwrite hazards)
                    float Ek = dst[dOff + k];
                    float Emk = dst[dOff + km];
                    float Ok = dst[dMid + k];
                    float Omk = dst[dMid + km];

                    // Correct Hartley pairwise mix:
                    // t (for index k):     c*Ok + s*Omk
                    // u (for index M - k): s*Ok - c*Omk  ← (fixed signs)
                    float t = c * Ok + s * Omk;
                    float u = s * Ok - c * Omk;

                    // Write four outputs
                    dst[dOff + k] = Ek + t;  // X[k]
                    dst[dOff + km] = Emk + u;  // X[M-k]
                    dst[dOff + M + k] = Ek - t;  // X[k+M]
                    dst[dOff + M + km] = Emk - u;  // X[M-k+M]
                }

                // Middle index k = M/2 if M is even (degenerate pair)
                if ((M & 1) == 0)
                {
                    int k = half;
                    float Ek = dst[dOff + k];
                    float Ok = dst[dMid + k];
                    dst[dOff + k] = Ek + Ok; // X[M/2]
                    dst[dOff + M + k] = Ek - Ok; // X[M + M/2]
                }
            }
            else
            {
                // Odd-length fallback: direct O(n^2) Hartley (pure, no FFT).
                DirectFHT(src, sOff, sStride, n, dst, dOff);
            }
        }

        /// <summary>
        /// Direct Hartley transform (no normalization): Y[k] = Σ_n x[n] * cas(2π n k / N).
        /// </summary>
        private static void DirectFHT(float[] src, int sOff, int sStride, int n, float[] dst, int dOff)
        {
            double w = 2.0 * Math.PI / n;

            for (int k = 0; k < n; k++)
            {
                double acc = 0.0;
                int p = sOff;
                for (int n0 = 0; n0 < n; n0++, p += sStride)
                {
                    acc += src[p] * Special.Cas((float)(w * n0 * k));
                }
                dst[dOff + k] = (float)acc;
            }
        }
        #endregion
    }
}
