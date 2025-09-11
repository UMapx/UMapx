using System;

namespace UMapx.Transform
{
    /// <summary>
    /// Fast real Weyl–Heisenberg (Gabor) transform via Hartley (FHT), matching the 2N×2N matrix (analysis exact).
    /// </summary>
    /// <remarks>
    /// Analysis factorization:
    ///   • Polyphase: n = a + b*M, N = M*L
    ///   • Correlation along L via Hartley-DFT (two FHT_L calls to recover C,S)
    ///   • For G2 (M/2 branch): for residues a ≥ M/2 apply extra +1 cyclic shift along L (implemented as spectral rotation)
    ///   • DFT along M via Hartley (two FHT_M calls → C,S) and quarter-period mixing by k mod 4
    /// Uses only your FastHartleyTransform; no complex arithmetic.
    /// Backward provided as exact O(N^2) matrix-form (Special.Cas).
    /// </remarks>
    [Serializable]
    public class FastRealWeylHeisenbergHartleyTransform : TransformBaseFloat, ITransform
    {
        #region Private data
        private readonly FastHartleyTransform FHT_L; // along L
        private readonly FastHartleyTransform FHT_M; // along M

        private readonly float[] g0;   // prototype window, length N
        private readonly int N, M, L;

        // quarter-period twiddles per k: cos(πk/2), sin(πk/2) ∈ {0,±1}
        private readonly sbyte[] c0;
        private readonly sbyte[] s0;

        // conj DFT of window polyphases along L, for each residue a
        private readonly float[][] Cg;   // Re spectrum for G_a
        private readonly float[][] Sg;   // Im-like component for G_a (sign matched for conj)
        private readonly float[][] Cg2;  // Re spectrum for G_{a+M/2}
        private readonly float[][] Sg2;  // Im-like component for G_{a+M/2}

        // cos/sin for rotation by e^{-i 2π q/L} (to apply +1 shift along L for a ≥ M/2 in G2)
        private readonly float[] cL;
        private readonly float[] sL;
        #endregion

        #region Ctor / Properties
        /// <summary>
        /// Fast real Weyl–Heisenberg (Gabor) transform via Hartley (FHT), matching the 2N×2N matrix (analysis exact).
        /// </summary>
        /// <param name="g0">Window function (orthogonalized)</param>
        /// <param name="M">Frequency shifts</param>
        /// <exception cref="ArgumentNullException">Exception</exception>
        /// <exception cref="ArgumentException">Exception</exception>
        public FastRealWeylHeisenbergHartleyTransform(float[] g0, int M)
        {
            if (g0 == null) throw new ArgumentNullException(nameof(g0));
            if (M <= 0 || (M & 1) != 0) throw new ArgumentException("M must be positive and even", nameof(M));

            this.N = g0.Length;
            if (N % M != 0) throw new ArgumentException("g0.Length must be divisible by M (N = M * L)", nameof(g0));

            this.M = M;
            this.L = N / M;
            this.g0 = (float[])g0.Clone();

            // Hartley transforms
            this.FHT_L = new FastHartleyTransform(false);
            this.FHT_M = new FastHartleyTransform(false);

            // quarter-period tables
            this.c0 = new sbyte[M];
            this.s0 = new sbyte[M];
            for (int k = 0; k < M; k++)
            {
                int r = k & 3;
                c0[k] = (sbyte)((r == 0) ? +1 : (r == 2) ? -1 : 0);
                s0[k] = (sbyte)((r == 1) ? +1 : (r == 3) ? -1 : 0);
            }

            // precompute cos/sin for e^{-i 2π q/L}
            this.cL = new float[L];
            this.sL = new float[L];
            for (int q = 0; q < L; q++)
            {
                double ang = 2.0 * Math.PI * q / L;
                cL[q] = (float)Math.Cos(ang);
                sL[q] = (float)Math.Sin(ang);
            }

            // window spectra along L
            this.Cg = new float[M][];
            this.Sg = new float[M][];
            this.Cg2 = new float[M][];
            this.Sg2 = new float[M][];

            PrecomputeWindowSpectra();
        }
        #endregion

        #region Private helpers

        /// <summary>
        /// Precompute “complex-like” Hartley spectra over L for all polyphase columns of the window:
        ///   G_a(q)  = Cg[a][q]  + i * Sg[a][q]
        ///   G2_a(q) = Cg2[a][q] + i * Sg2[a][q], where a' = a + M/2 (mod M).
        ///
        /// We obtain Re/Im-like parts via two FHT_L calls (vector and its cyclic-reverse):
        ///   Re = 0.5 * (H + H_rev),  Im = 0.5 * (H - H_rev)
        ///
        /// Notes:
        /// • The sign convention for S(·) is chosen so that multiplying by conj(G) in the
        ///   Hartley domain maps to the (re − im) combination later in Forward.
        /// • FHT normalization follows your FastHartleyTransform settings (normalized or not).
        /// </summary>
        private void PrecomputeWindowSpectra()
        {
            // Polyphase buffers (length L)
            var Ga = new float[L];
            var G2a = new float[L];

            for (int a = 0; a < M; a++)
            {
                // -----------------------------
                // residue a → spectra Cg[a], Sg[a]
                // -----------------------------
                // Build polyphase column: Ga[b] = g0[a + b*M], b = 0..L-1
                for (int b = 0; b < L; b++) Ga[b] = g0[a + b * M];

                // One FHT_L; the spectrum of the cyclic-reversed signal is H[(L-q) % L]
                var H = FHT_L.Forward(Ga);

                var C = new float[L];
                var S = new float[L];

                // Recover cos/sin sums:
                //   C[q] = 0.5 * (H[q] + H_rev[q])  = Σ g * cos
                //   S[q] = 0.5 * (H[q] - H_rev[q])  = Σ g * sin
                // where H_rev[q] = H[(L - q) % L]
                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);  // mirrors 0→0, (L/2)→(L/2) when L is even
                    float HR = H[qrev];
                    C[q] = 0.5f * (H[q] + HR); // Re
                    S[q] = 0.5f * (H[q] - HR); // Im-like (sign matched for conj)
                }

                Cg[a] = C;
                Sg[a] = S;

                // -----------------------------
                // residue a' = a + M/2 → spectra Cg2[a], Sg2[a]
                // -----------------------------
                int ap = (a + (M >> 1)) % M;

                for (int b = 0; b < L; b++) G2a[b] = g0[ap + b * M];

                var H2 = FHT_L.Forward(G2a);

                var C2 = new float[L];
                var S2 = new float[L];

                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);
                    float HR = H2[qrev];
                    C2[q] = 0.5f * (H2[q] + HR);
                    S2[q] = 0.5f * (H2[q] - HR);
                }

                Cg2[a] = C2;
                Sg2[a] = S2;
            }
        }

        #endregion

        #region Transform methods

        /// <summary>
        /// Forward (analysis): y[2N] -> c[2N] equals G^T * y (exactly matches your 2N×2N matrix).
        /// If y.Length == N, the bottom half (xi) is assumed to be zeros.
        /// </summary>
        /// <remarks>
        /// Factorization:
        /// 1) Along L (time steps): compute correlations r_a[l], s_a[l] with window polyphases,
        ///    entirely in the Hartley domain using two FHT_L calls per vector to extract cos/sin.
        ///    For branch G2, for residues a ≥ M/2, apply an extra +1 shift along L, implemented
        ///    as spectral rotation by e^{-i 2π q / L}.
        /// 2) Along M (frequency shifts): FHT_M + reverse trick to get cos/sin for each of r,s,
        ///    then quarter-period mixing using k mod 4 via c0[k]=cos(πk/2), s0[k]=sin(πk/2) ∈ {0,±1}.
        /// 3) Branch formulas:
        ///      c1 =  cos(φ − πk/2)·R  +  sin(φ − πk/2)·S
        ///      c2 = −sin(φ − πk/2)·R2 +  cos(φ − πk/2)·S2
        /// 4) A final scale = 1/L is applied to match your matrix convention.
        /// </remarks>
        public override float[] Forward(float[] y)
        {
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (y.Length != N && y.Length != 2 * N)
                throw new ArgumentException($"Input length must be N={N} or 2N={2 * N}");

            // split input into top/bottom halves (xr, xi)
            var xr = new float[N];
            var xi = new float[N];

            if (y.Length == 2 * N)
            {
                Buffer.BlockCopy(y, 0, xr, 0, sizeof(float) * N);
                Buffer.BlockCopy(y, sizeof(float) * N, xi, 0, sizeof(float) * N);
            }
            else
            {
                Buffer.BlockCopy(y, 0, xr, 0, sizeof(float) * N);
                // xi stays zero
            }

            // 1) Correlation along L for each residue a (for both branches and both halves)
            var R = new float[M][];
            var S = new float[M][];
            var R2 = new float[M][];
            var S2 = new float[M][];

            var Xa = new float[L];
            var Hy = new float[L];
            var Hy2 = new float[L];

            for (int a = 0; a < M; a++)
            {
                // --- xr with Ga ---
                for (int b = 0; b < L; b++) Xa[b] = xr[a + b * M];

                // Single FHT_L; cos/sin via spectral mirroring: HR[q] = H[(L-q)%L]
                var Hx = FHT_L.Forward(Xa);

                var Cx = new float[L];
                var Sx = new float[L];
                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);
                    float HR = Hx[qrev];
                    Cx[q] = 0.5f * (Hx[q] + HR); // Σ xr cos
                    Sx[q] = 0.5f * (Hx[q] - HR); // Σ xr sin
                }

                // --- xi with Ga ---
                for (int b = 0; b < L; b++) Xa[b] = xi[a + b * M];

                var Hi = FHT_L.Forward(Xa);

                var Cxi = new float[L];
                var Sxi = new float[L];
                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);
                    float HR = Hi[qrev];
                    Cxi[q] = 0.5f * (Hi[q] + HR); // Σ xi cos
                    Sxi[q] = 0.5f * (Hi[q] - HR); // Σ xi sin
                }

                // Y = X * conj(G) for Ga
                var Cg_a = Cg[a]; var Sg_a = Sg[a];

                for (int q = 0; q < L; q++)
                {
                    float re = Cx[q] * Cg_a[q] + Sx[q] * Sg_a[q];
                    float im = Cx[q] * Sg_a[q] - Sx[q] * Cg_a[q];
                    Hy[q] = re - im;   // xr with Ga

                    float re_i = Cxi[q] * Cg_a[q] + Sxi[q] * Sg_a[q];
                    float im_i = Cxi[q] * Sg_a[q] - Sxi[q] * Cg_a[q];
                    Hy2[q] = re_i - im_i; // xi with Ga
                }

                R[a] = FHT_L.Backward(Hy);   // r_a[l]
                S[a] = FHT_L.Backward(Hy2);  // s_a[l]

                // Y = X * conj(G2) for Ga2 (residue a+M/2), with extra +1 shift along L when a ≥ M/2
                var Cg2_a = Cg2[a]; var Sg2_a = Sg2[a];
                bool shift1 = a >= (M >> 1);

                for (int q = 0; q < L; q++)
                {
                    // rotate conj(G2) by e^{-i·2π q/L} if shift1
                    float Cr = Cg2_a[q];
                    float Sr = Sg2_a[q];

                    if (shift1)
                    {
                        float cx = cL[q], s = sL[q];
                        float Cr2 = Cr * cx + Sr * s;
                        float Sr2 = -Cr * s + Sr * cx;
                        Cr = Cr2; Sr = Sr2;
                    }

                    float re2 = Cx[q] * Cr + Sx[q] * Sr;
                    float im2 = Cx[q] * Sr - Sx[q] * Cr;
                    Hy[q] = re2 - im2;

                    float re2i = Cxi[q] * Cr + Sxi[q] * Sr;
                    float im2i = Cxi[q] * Sr - Sxi[q] * Cr;
                    Hy2[q] = re2i - im2i;
                }

                R2[a] = FHT_L.Backward(Hy);   // r2_a[l]
                S2[a] = FHT_L.Backward(Hy2);  // s2_a[l]
            }

            // 2) For each l: FHT along a (length M) + quarter-period mixing
            var c = new float[2 * N];

            var vec = new float[M];

            var CR = new float[M]; var SR = new float[M];
            var CS = new float[M]; var SS = new float[M];

            // scale for M-stage mixing:
            float scale = 1.0f / L;

            for (int l = 0; l < L; l++)
            {
                // --- branch 1 (G1): from R,S ---
                for (int a = 0; a < M; a++) vec[a] = R[a][l];
                var H1 = FHT_M.Forward(vec);
                for (int k = 0; k < M; k++)
                {
                    int krev = (k == 0) ? 0 : (M - k);
                    float HR = H1[krev];
                    CR[k] = 0.5f * (H1[k] + HR);
                    SR[k] = 0.5f * (H1[k] - HR);
                }

                for (int a = 0; a < M; a++) vec[a] = S[a][l];
                H1 = FHT_M.Forward(vec);
                for (int k = 0; k < M; k++)
                {
                    int krev = (k == 0) ? 0 : (M - k);
                    float HR = H1[krev];
                    CS[k] = 0.5f * (H1[k] + HR);
                    SS[k] = 0.5f * (H1[k] - HR);
                }

                for (int k = 0; k < M; k++)
                {
                    // c1 = cos(φ−πk/2)·R + sin(φ−πk/2)·S
                    float re1 = c0[k] * (CR[k] + SS[k]) + s0[k] * (SR[k] - CS[k]);
                    c[l * M + k] = re1 * scale;
                }

                // --- branch 2 (G2): from R2,S2) ---
                for (int a = 0; a < M; a++) vec[a] = R2[a][l];
                H1 = FHT_M.Forward(vec);
                for (int k = 0; k < M; k++)
                {
                    int krev = (k == 0) ? 0 : (M - k);
                    float HR = H1[krev];
                    CR[k] = 0.5f * (H1[k] + HR);
                    SR[k] = 0.5f * (H1[k] - HR);
                }

                for (int a = 0; a < M; a++) vec[a] = S2[a][l];
                H1 = FHT_M.Forward(vec);
                for (int k = 0; k < M; k++)
                {
                    int krev = (k == 0) ? 0 : (M - k);
                    float HR = H1[krev];
                    CS[k] = 0.5f * (H1[k] + HR);
                    SS[k] = 0.5f * (H1[k] - HR);
                }

                for (int k = 0; k < M; k++)
                {
                    // c2 = −sin(φ−πk/2)·R2 + cos(φ−πk/2)·S2
                    float re2 = c0[k] * (-SR[k] + CS[k]) + s0[k] * (CR[k] + SS[k]);
                    c[l * M + k + N] = re2 * scale;
                }
            }

            return c;
        }

        /// <summary>
        /// Backward (synthesis): c[2N] -> y[2N], fast Hartley-only inverse
        /// that mirrors the Forward factorization (matches your matrix form).
        /// </summary>
        /// <remarks>
        /// Inverse steps:
        /// 1) For each l, recover four k-sums via 2×FHT_M with quarter-period tables:
        ///    A1 = Σ c1·cos(φ−πk/2),  A2 = Σ c1·sin(φ−πk/2),
        ///    B1 = Σ c2·cos(φ−πk/2),  B2 = Σ c2·sin(φ−πk/2).
        ///    We do **one** FHT per weighted vector and obtain the “reverse” via spectral
        ///    mirroring (H_rev[k] = H[(M-k) mod M]).
        /// 2) For each residue a, do two L-convolutions in the Hartley domain:
        ///    xr: y1 = g_a * A1 ; y2 = g_{a+M/2}(+1 shift if a≥M/2) * (−B2)
        ///    xi: z1 = g_a * A2 ; z2 = g_{a+M/2}(+1 shift if a≥M/2) * (+B1)
        ///    Each convolution uses **one** forward FHT_L (cos/sin via mirroring) and
        ///    one inverse FHT_L.
        /// 3) Inverse FHT_L and deinterleave to n = a + b*M into xr and xi.
        /// </remarks>
        public override float[] Backward(float[] c)
        {
            if (c == null) throw new ArgumentNullException(nameof(c));
            if (c.Length != 2 * N) throw new ArgumentException($"Input length must be 2N = {2 * N}");

            var xr = new float[N];
            var xi = new float[N];

            // Undo the Forward mixing scale (Forward stored re*scale, here we multiply by 1/scale = L).
            float invScaleM = 1.0f / L;

            // === Step 1: recover A1, A2, B1, B2 for each l via two FHT_M calls ===
            var A1 = new float[M][]; var A2 = new float[M][];
            var B1 = new float[M][]; var B2 = new float[M][];
            for (int a = 0; a < M; a++) { A1[a] = new float[L]; A2[a] = new float[L]; B1[a] = new float[L]; B2[a] = new float[L]; }

            var C1 = new float[M]; var C2 = new float[M]; // c1(l,·), c2(l,·) with scale undone
            var P = new float[M]; var Q = new float[M]; // c0-weighted and s0-weighted

            for (int l = 0; l < L; l++)
            {
                // Unpack c1, c2 for fixed l and remove the Forward scale
                for (int k = 0; k < M; k++)
                {
                    int u = l * M + k;
                    C1[k] = c[u] * invScaleM;
                    C2[k] = c[u + N] * invScaleM;
                }

                // --- A1 = cos_sum(c0*C1) + sin_sum(s0*C1) ---
                for (int k = 0; k < M; k++) { P[k] = c0[k] * C1[k]; Q[k] = s0[k] * C1[k]; }

                var HP = FHT_M.Forward(P);
                var HQ = FHT_M.Forward(Q);

                for (int a = 0; a < M; a++)
                {
                    int arev = (a == 0) ? 0 : (M - a);
                    float cosP = 0.5f * (HP[a] + HP[arev]); // cos_sum(P)
                    float sinQ = 0.5f * (HQ[a] - HQ[arev]); // sin_sum(Q)
                    A1[a][l] = cosP + sinQ;
                }

                // --- A2 = sin_sum(c0*C1) − cos_sum(s0*C1) ---
                for (int a = 0; a < M; a++)
                {
                    int arev = (a == 0) ? 0 : (M - a);
                    float cosS = 0.5f * (HQ[a] + HQ[arev]); // cos_sum(s0*C1)
                    float sinC = 0.5f * (HP[a] - HP[arev]); // sin_sum(c0*C1)
                    A2[a][l] = sinC - cosS;
                }

                // --- B1 = cos_sum(c0*C2) + sin_sum(s0*C2) ---
                for (int k = 0; k < M; k++) { P[k] = c0[k] * C2[k]; Q[k] = s0[k] * C2[k]; }

                HP = FHT_M.Forward(P);
                HQ = FHT_M.Forward(Q);

                for (int a = 0; a < M; a++)
                {
                    int arev = (a == 0) ? 0 : (M - a);
                    float cosP = 0.5f * (HP[a] + HP[arev]); // cos_sum(P)
                    float sinQ = 0.5f * (HQ[a] - HQ[arev]); // sin_sum(Q)
                    B1[a][l] = cosP + sinQ;
                }

                // --- B2 = sin_sum(c0*C2) − cos_sum(s0*C2) ---
                for (int a = 0; a < M; a++)
                {
                    int arev = (a == 0) ? 0 : (M - a);
                    float cosS_fix = 0.5f * (HQ[a] + HQ[arev]); // cos_sum(s0*C2)
                    float sinC_fix = 0.5f * (HP[a] - HP[arev]); // sin_sum(c0*C2)
                    B2[a][l] = sinC_fix - cosS_fix;
                }
            }

            // === Step 2: two L-convolutions per residue a, then deinterleave to n = a + b*M ===
            var col = new float[L];
            var Cx = new float[L];
            var Sx = new float[L];
            var Hy = new float[L];

            for (int a = 0; a < M; a++)
            {
                // ----- xr: y1 = g_a * A1 -----
                Array.Copy(A1[a], col, L);
                var Hx = FHT_L.Forward(col);
                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);
                    float HR = Hx[qrev];
                    Cx[q] = 0.5f * (Hx[q] + HR);
                    Sx[q] = 0.5f * (Hx[q] - HR);
                }

                var Cg_a = Cg[a]; var Sg_a = Sg[a];
                for (int q = 0; q < L; q++)
                {
                    float re = Cx[q] * Cg_a[q] - Sx[q] * Sg_a[q];
                    float im = Cx[q] * Sg_a[q] + Sx[q] * Cg_a[q];
                    Hy[q] = re + im; // Hartley spectrum (cas)
                }
                var y1 = FHT_L.Backward(Hy);

                // ----- xr: y2 = g_{a+M/2}(+1 shift if a≥M/2) * (−B2) -----
                Array.Copy(B2[a], col, L);
                for (int i = 0; i < L; i++) col[i] = -col[i]; // minus sign from the branch-2 formula
                Hx = FHT_L.Forward(col);
                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);
                    float HR = Hx[qrev];
                    Cx[q] = 0.5f * (Hx[q] + HR);
                    Sx[q] = 0.5f * (Hx[q] - HR);
                }

                var Cg2_a = Cg2[a]; var Sg2_a = Sg2[a];
                bool shift1 = a >= (M >> 1);
                for (int q = 0; q < L; q++)
                {
                    float Cr = Cg2_a[q], Sr = Sg2_a[q];
                    if (shift1)
                    {
                        // +1 shift along L => multiply G2 by e^{-i 2π q / L}
                        float cph = cL[q], sph = sL[q];
                        float Cr2 = Cr * cph + Sr * sph;
                        float Sr2 = -Cr * sph + Sr * cph;
                        Cr = Cr2; Sr = Sr2;
                    }
                    float re = Cx[q] * Cr - Sx[q] * Sr;
                    float im = Cx[q] * Sr + Sx[q] * Cr;
                    Hy[q] = re + im;
                }
                var y2 = FHT_L.Backward(Hy);

                // Scatter back to xr at n = a + b*M
                for (int b = 0; b < L; b++)
                {
                    int n = a + b * M;
                    xr[n] = y1[b] + y2[b];
                }

                // ----- xi: z1 = g_a * A2 -----
                Array.Copy(A2[a], col, L);
                Hx = FHT_L.Forward(col);
                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);
                    float HR = Hx[qrev];
                    Cx[q] = 0.5f * (Hx[q] + HR);
                    Sx[q] = 0.5f * (Hx[q] - HR);
                }
                for (int q = 0; q < L; q++)
                {
                    float re = Cx[q] * Cg_a[q] - Sx[q] * Sg_a[q];
                    float im = Cx[q] * Sg_a[q] + Sx[q] * Cg_a[q];
                    Hy[q] = re + im;
                }
                var z1 = FHT_L.Backward(Hy);

                // ----- xi: z2 = g_{a+M/2}(+1 shift if a≥M/2) * (+B1) -----
                Array.Copy(B1[a], col, L);
                Hx = FHT_L.Forward(col);
                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);
                    float HR = Hx[qrev];
                    Cx[q] = 0.5f * (Hx[q] + HR);
                    Sx[q] = 0.5f * (Hx[q] - HR);
                }
                for (int q = 0; q < L; q++)
                {
                    float Cr = Cg2_a[q], Sr = Sg2_a[q];
                    if (shift1)
                    {
                        // same e^{-i 2π q / L} rotation
                        float cph = cL[q], sph = sL[q];
                        float Cr2 = Cr * cph + Sr * sph;
                        float Sr2 = -Cr * sph + Sr * cph;
                        Cr = Cr2; Sr = Sr2;
                    }
                    float re = Cx[q] * Cr - Sx[q] * Sr;
                    float im = Cx[q] * Sr + Sx[q] * Cr;
                    Hy[q] = re + im;
                }
                var z2 = FHT_L.Backward(Hy);

                for (int b = 0; b < L; b++)
                {
                    int n = a + b * M;
                    xi[n] = z1[b] + z2[b];
                }
            }

            // Pack y = [xr; xi]
            var y = new float[2 * N];
            Buffer.BlockCopy(xr, 0, y, 0, sizeof(float) * N);
            Buffer.BlockCopy(xi, 0, y, sizeof(float) * N, sizeof(float) * N);
            return y;
        }

        #endregion
    }
}
