using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Hartley transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
    /// </remarks>
    [Serializable]
    public class FastHartleyTransform : TransformBaseFloat, ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FastFourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Hartley transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="spectrumType">Spectrum type</param>
        public FastHartleyTransform(bool normalized = true, Direction direction = Direction.Vertical, SpectrumType spectrumType = SpectrumType.Fourier)
        {
            this.FFT = new FastFourierTransform(normalized, direction);
            this.SpectrumType = spectrumType;
        }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        /// <summary>
        /// Gets or sets the processing direction.
        /// </summary>
        public override Direction Direction
        {
            get
            {
                return this.FFT.Direction;
            }
            set
            {
                this.FFT.Direction = value;
            }
        }
        /// <summary>
        /// Gets or sets spectrum type.
        /// </summary>
        public SpectrumType SpectrumType { get; set; }
        #endregion

        #region Fast Hartley Transform (1D, float)
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override float[] Forward(float[] A)
        {
            if (SpectrumType == SpectrumType.Fourier)
            {
                Complex32[] B = Matrice.ToComplex(A);
                B = FFT.Forward(B);

                int length = A.Length, i;
                float[] Hk = new float[length];

                for (i = 0; i < length; i++)
                {
                    Hk[i] = B[i].Real - B[i].Imag;
                }

                return Hk;
            }
            else
            {
                int N = A.Length;
                var outv = new float[N];
                FHT(A, 0, 1, N, outv, 0);

                if (Normalized)
                {
                    float s = 1f / Maths.Sqrt(N);
                    for (int i = 0; i < N; i++) outv[i] *= s;
                }
                return outv;
            }
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override float[] Backward(float[] B)
        {
            if (SpectrumType == SpectrumType.Fourier)
            {
                Complex32[] A = Matrice.ToComplex(B);
                A = FFT.Backward(A);

                int length = B.Length, i;
                float[] Hk = new float[length];

                for (i = 0; i < length; i++)
                {
                    Hk[i] = A[i].Real + A[i].Imag;
                }

                return Hk;
            }
            else
            {
                int N = B.Length;
                var outv = new float[N];
                FHT(B, 0, 1, N, outv, 0);

                if (Normalized)
                {
                    float s = 1f / Maths.Sqrt(N);
                    for (int i = 0; i < N; i++) outv[i] *= s;
                }
                return outv;
            }
        }
        #endregion

        #region Core FHT
        /// <summary>
        /// Radix-2 DIT FHT for even n; safe O(n^2) direct Hartley for odd n.
        /// Writes n outputs contiguously to dst[dOff .. dOff + n - 1].
        /// Optimizations:
        ///  • Fast bases for n=2 and n=4 (no recursion, no trig)
        ///  • Trigonometric recurrences for twiddles (one sin/cos per call, O(1) per k)
        ///  • Single-precision math (MathF) to reduce casts and overhead
        /// </summary>
        private static void FHT(float[] src, int sOff, int sStride, int n, float[] dst, int dOff)
        {
            if (n <= 0) return;

            // Base cases (cut recursion + trig calls):
            if (n == 1)
            {
                dst[dOff] = src[sOff];
                return;
            }
            if (n == 2)
            {
                float a = src[sOff];
                float b = src[sOff + sStride];
                dst[dOff + 0] = a + b;
                dst[dOff + 1] = a - b;
                return;
            }
            if (n == 4)
            {
                int s1 = sOff + sStride;
                int s2 = s1 + sStride;
                int s3 = s2 + sStride;

                float a = src[sOff];
                float b = src[s1];
                float c = src[s2];
                float d = src[s3];

                // 4-point DHT (Hartley) matrix:
                // [1  1  1  1;
                //  1  1 -1 -1;
                //  1 -1  1 -1;
                //  1 -1 -1  1]
                dst[dOff + 0] = (a + b) + (c + d);
                dst[dOff + 1] = (a + b) - (c + d);
                dst[dOff + 2] = (a - b) + (c - d);
                dst[dOff + 3] = (a - b) - (c - d);
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

                // k = 0 (degenerate pair): no trig here
                {
                    float E0 = dst[dOff + 0];
                    float O0 = dst[dMid + 0];
                    dst[dOff + 0] = E0 + O0; // X[0]
                    dst[dOff + M + 0] = E0 - O0; // X[M]
                }

                // Pairwise combine with trig recurrences:
                // We need c_k = cos(πk/M), s_k = sin(πk/M), k=1..half-1.
                int half = M >> 1;

                if (half > 1)
                {
                    // One sin/cos for the step θ = π/M
                    float theta = Maths.Pi / M;
                    float cs = Maths.Cos(theta);
                    float ss = Maths.Sin(theta);

                    // Start at k = 1: (c, s) = (cos θ, sin θ)
                    float c = cs;
                    float s = ss;

                    for (int k = 1; k < half; k++)
                    {
                        int km = M - k;

                        // Read inputs (avoid overwrite hazards)
                        float Ek = dst[dOff + k];
                        float Emk = dst[dOff + km];
                        float Ok = dst[dMid + k];
                        float Omk = dst[dMid + km];

                        // Canonical Hartley pairwise mix:
                        // t (for index  k)  =  c*Ok + s*Omk
                        // u (for index M-k) =  s*Ok - c*Omk
                        float t = c * Ok + s * Omk;
                        float u = s * Ok - c * Omk;

                        dst[dOff + k] = Ek + t;  // X[k]
                        dst[dOff + km] = Emk + u;  // X[M-k]
                        dst[dOff + M + k] = Ek - t;  // X[k+M]
                        dst[dOff + M + km] = Emk - u;  // X[M-k+M]

                        // Update (c, s) → (cos((k+1)θ), sin((k+1)θ)) via rotation:
                        // c' = c*cs - s*ss;  s' = s*cs + c*ss
                        float cNext = c * cs - s * ss;
                        s = s * cs + c * ss;
                        c = cNext;
                    }
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
        /// Optimized with trigonometric recurrences per k (no per-sample trig calls).
        /// </summary>
        private static void DirectFHT(float[] src, int sOff, int sStride, int n, float[] dst, int dOff)
        {
            // For each k, cas(m * ωk) is generated via recurrence:
            // cos((m+1)φ) = cos mφ * cos φ - sin mφ * sin φ
            // sin((m+1)φ) = sin mφ * cos φ + cos mφ * sin φ
            // then cas = cos + sin.
            float twoPiOverN = 2f * Maths.Pi / n;

            for (int k = 0; k < n; k++)
            {
                float phi = twoPiOverN * k;       // φ = 2π k / N
                float cStep = Maths.Cos(phi);     // cos φ
                float sStep = Maths.Sin(phi);     // sin φ

                // m = 0 → cos 0 = 1, sin 0 = 0
                float c = 1f;
                float s = 0f;

                // params
                float acc = 0f;
                int p = sOff;

                for (int m = 0; m < n; m++, p += sStride)
                {
                    // cas(mφ) = cos(mφ) + sin(mφ) = c + s
                    acc += src[p] * (c + s);

                    // Advance (c, s) to (m+1) using one rotation
                    // c' = c*cStep - s*sStep;   s' = s*cStep + c*sStep
                    float cNext = c * cStep - s * sStep;
                    s = s * cStep + c * sStep;
                    c = cNext;
                }

                dst[dOff + k] = acc;
            }
        }
        #endregion

    }
}
