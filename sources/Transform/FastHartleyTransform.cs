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
        /// <param name="spectrumType">Spectrum type</param>
        /// <param name="direction">Processing direction</param>
        public FastHartleyTransform(bool normalized = true, SpectrumType spectrumType = SpectrumType.Fourier, Direction direction = Direction.Vertical)
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

        #region Fast Hartley Transform
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

                // 2-point DHT (Hartley) matrix:
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
                FHT_BluesteinHartley(src, sOff, sStride, n, dst, dOff);
            }
        }

        /// <summary>
        /// Pure-Hartley Bluestein (no FFT):
        /// DHT_N(x)[k] = Re{X[k]} - Im{X[k]}, where
        /// X[k] = e^{-jπ k^2/N} · ( a * b )[k],
        /// a[n] = x[n] · e^{-jπ n^2/N},   b[m] = e^{+jπ m^2/N}.
        /// The linear convolution (a*b) is computed via FHT-only engine:
        /// for real sequences, conv = IDHT( combine(Ha, Hb) ), where
        /// combine implements the DHT convolution theorem; complex conv
        /// is obtained by four real convs.
        /// </summary>
        private static void FHT_BluesteinHartley(float[] src, int sOff, int sStride, int n, float[] dst, int dOff)
        {
            // 1) Build chirped a = x * exp(-j π n^2 / N)
            var ar = new float[n];
            var ai = new float[n];

            float c0 = Maths.Pi / n; // π / N
            int p = sOff;
            for (int i = 0; i < n; i++, p += sStride)
            {
                float x = src[p];
                float ang = c0 * i * i;      // π i^2 / N
                float c = Maths.Cos(ang);    // cos(ang)
                float s = -Maths.Sin(ang);   // -sin(ang) for exp(-j ang)

                ar[i] = x * c;
                ai[i] = x * s;
            }

            // 2) Build b (kernel) for linear convolution on length P ≥ 2N-1.
            int L = 2 * n - 1;
            int P = NextPow2(L);

            var br = new float[P]; // real part of b
            var bi = new float[P]; // imag part of b

            // Positive indices m = 0..N-1
            for (int m = 0; m < n; m++)
            {
                float ang = c0 * m * m;  // π m^2 / N
                br[m] = Maths.Cos(ang);
                bi[m] = Maths.Sin(ang);  // +sin for exp(+j ang)
            }
            // Negative indices m = -(N-1)..-1 mapped to P-(N-1)..P-1
            for (int mm = 1; mm <= n - 1; mm++)
            {
                int idx = P - mm;        // corresponds to m = -mm
                float mneg = mm;         // (-mm)^2 = mm^2
                float ang = c0 * mneg * mneg;
                br[idx] = Maths.Cos(ang);
                bi[idx] = Maths.Sin(ang);
            }

            // 3) Compute complex linear convolution c = a (*) b via Hartley-only engine
            //    Using: (ar+jai)*(br+jbi) = (ar⊗br - ai⊗bi) + j(ar⊗bi + ai⊗br),
            //    where ⊗ is linear convolution computed via FHT on size P.

            // Real convs (length L, calculated on size P):
            var t1 = RealLinearConvFHT(ar, n, br, P, L); // ar ⊗ br
            var t2 = RealLinearConvFHT(ai, n, bi, P, L); // ai ⊗ bi
            var t3 = RealLinearConvFHT(ar, n, bi, P, L); // ar ⊗ bi
            var t4 = RealLinearConvFHT(ai, n, br, P, L); // ai ⊗ br

            var cr = new float[L];
            var ci = new float[L];
            for (int k = 0; k < L; k++)
            {
                cr[k] = t1[k] - t2[k];
                ci[k] = t3[k] + t4[k];
            }

            // 4) X[k] = e^{-jπ k^2/N} * c[k], take k=0..N-1
            //    Then DHT = Re(X) - Im(X).
            for (int k = 0; k < n; k++)
            {
                float ang = c0 * k * k;  // π k^2 / N
                float c = Maths.Cos(ang);
                float s = -Maths.Sin(ang); // for exp(-j ang)

                float rr = cr[k];
                float ii = ci[k];

                // Complex multiply (rr + j*ii) * (c + j*s)
                float xr = rr * c - ii * s;
                float xi = rr * s + ii * c;

                dst[dOff + k] = xr - xi; // Hartley: Re - Im
            }
        }

        /// <summary>
        /// Real linear convolution y = a ⊗ b (length L = la + lb_eff - 1),
        /// computed via Hartley-only convolution theorem on size P (power-of-two).
        /// Here b is already padded on size P (lb_eff = P), and a occupies 0..la-1.
        /// Implementation details:
        ///  1) Zero-pad a and b to length P
        ///  2) Y_H = HartleyConvolutionCombine( FHT(a), FHT(b) )         (circular on P)
        ///  3) y_pad = IDHT(Y_H) / P                                     (time domain)
        ///  4) Return first L samples (linear conv region)
        /// Note: uses pure FHT (no normalization inside; 1/P appears explicitly).
        /// </summary>
        private static float[] RealLinearConvFHT(float[] a, int la, float[] b_pad, int P, int L)
        {
            // Safety: Pchecked==P is just to stress size is power-of-two
            // Pack a→A (size P), b already in b_pad (size P)
            var A = new float[P];
            var B = new float[P];
            for (int i = 0; i < la; i++) A[i] = a[i];
            // b_pad may carry data at both head and tail to emulate negative indices
            Array.Copy(b_pad, 0, B, 0, P);

            // FHT(A), FHT(B) on size P (pure Hartley, radix-2)
            var AH = new float[P];
            var BH = new float[P];
            FHT(A, 0, 1, P, AH, 0);
            FHT(B, 0, 1, P, BH, 0);

            // Combine in Hartley domain (circular conv theorem for DHT):
            // Let rev(k) = (P - k) % P. If F = H(a), G = H(b),
            // Then H(a ⊛ b)[k] = 1/2 * ( F[k]*G[k] + F[rev]*G[k] + F[k]*G[rev] - F[rev]*G[rev] )
            // Proof: via mapping DFT↔DHT (Re/Im vs H/H_rev), per-k pair algebra.
            var YH = new float[P];
            for (int k = 0; k < P; k++)
            {
                int rk = (k == 0) ? 0 : (P - k);
                float Fk = AH[k];
                float Fr = AH[rk];
                float Gk = BH[k];
                float Gr = BH[rk];

                YH[k] = 0.5f * (Fk * Gk + Fr * Gk + Fk * Gr - Fr * Gr);
            }

            // y_pad = IDHT(YH) = FHT(YH) / P  (since unnormalized FHT is self-inverse up to factor P)
            var ypad = new float[P];
            FHT(YH, 0, 1, P, ypad, 0);

            float invP = 1f / P;
            for (int i = 0; i < P; i++) ypad[i] *= invP;

            // Return first L samples (valid linear convolution region)
            var y = new float[L];
            Array.Copy(ypad, 0, y, 0, L);
            return y;
        }

        /// <summary>
        /// Next power of two ≥ n (n ≥ 1).
        /// </summary>
        private static int NextPow2(int n)
        {
            int p = 1;
            while (p < n) p <<= 1;
            return p;
        }
        #endregion
    }
}
