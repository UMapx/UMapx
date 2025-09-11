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
        /// Hartley-specific “cyclic reverse” used to separate cos/sin sums from a single FHT.
        /// It swaps v[1] ↔ v[n-1], v[2] ↔ v[n-2], … while keeping v[0] intact.
        /// With this permutation:
        ///   H_rev = FHT( v_rev )
        /// and therefore:
        ///   cos_sum = 0.5 * (H + H_rev)
        ///   sin_sum = 0.5 * (H - H_rev)
        /// </summary>
        private static void CyclicReverseInPlace(float[] v)
        {
            int n = v.Length, half = n / 2;
            for (int i = 1; i <= half; i++)
            {
                int a = i, b = n - i;
                var t = v[a]; v[a] = v[b]; v[b] = t;
            }
        }

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
            var Ga = new float[L];
            var G2a = new float[L];
            var tmp = new float[L];

            for (int a = 0; a < M; a++)
            {
                // Polyphase column for residue a: g0[a + b*M], b = 0..L-1
                for (int b = 0; b < L; b++) Ga[b] = g0[a + b * M];

                // Two FHTs to split into cos/sin
                var H = FHT_L.Forward(Ga);
                Array.Copy(Ga, tmp, L); CyclicReverseInPlace(tmp);
                var HR = FHT_L.Forward(tmp);

                var C = new float[L];
                var S = new float[L];
                for (int q = 0; q < L; q++)
                {
                    C[q] = 0.5f * (H[q] + HR[q]); // Re-like
                    S[q] = 0.5f * (H[q] - HR[q]); // Im-like (sign aligned to conj usage)
                }
                Cg[a] = C; Sg[a] = S;

                // Same for residue a' = a + M/2 (mod M)
                int ap = (a + (M >> 1)) % M;
                for (int b = 0; b < L; b++) G2a[b] = g0[ap + b * M];

                var H2 = FHT_L.Forward(G2a);
                Array.Copy(G2a, tmp, L); CyclicReverseInPlace(tmp);
                var H2R = FHT_L.Forward(tmp);

                var C2 = new float[L];
                var S2 = new float[L];
                for (int q = 0; q < L; q++)
                {
                    C2[q] = 0.5f * (H2[q] + H2R[q]);
                    S2[q] = 0.5f * (H2[q] - H2R[q]);
                }
                Cg2[a] = C2; Sg2[a] = S2;
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

            // Split input into top/bottom halves (xr, xi), each of length N
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
                // xi remains zero
            }

            // 1) Correlation along L for each residue a (both branches, both halves)
            var R = new float[M][];
            var S = new float[M][];
            var R2 = new float[M][];
            var S2 = new float[M][];

            var Xa = new float[L]; // polyphase slice of xr for residue a
            var Xrv = new float[L]; // its cyclic-reverse (to extract sin-sum)
            var Hy = new float[L]; // Hartley spectrum of correlation for xr
            var Hy2 = new float[L]; // Hartley spectrum of correlation for xi

            for (int a = 0; a < M; a++)
            {
                // --- xr with G_a ---
                for (int b = 0; b < L; b++) Xa[b] = xr[a + b * M];

                // Two FHTs to get cos/sin sums along L
                var Hx = FHT_L.Forward(Xa);
                Array.Copy(Xa, Xrv, L); CyclicReverseInPlace(Xrv);
                var HxR = FHT_L.Forward(Xrv);

                var Cx = new float[L];
                var Sx = new float[L];
                for (int q = 0; q < L; q++)
                {
                    Cx[q] = 0.5f * (Hx[q] + HxR[q]); // Σ xr cos
                    Sx[q] = 0.5f * (Hx[q] - HxR[q]); // Σ xr sin
                }

                // --- xi with G_a ---
                for (int b = 0; b < L; b++) Xa[b] = xi[a + b * M];

                var Hi = FHT_L.Forward(Xa);
                Array.Copy(Xa, Xrv, L); CyclicReverseInPlace(Xrv);
                var HiR = FHT_L.Forward(Xrv);

                var Cxi = new float[L];
                var Sxi = new float[L];
                for (int q = 0; q < L; q++)
                {
                    Cxi[q] = 0.5f * (Hi[q] + HiR[q]); // Σ xi cos
                    Sxi[q] = 0.5f * (Hi[q] - HiR[q]); // Σ xi sin
                }

                // Multiply by conj(G_a) in Hartley domain:
                // (Cx + iSx) * (Cg − iSg) → Hy = (re − im) combination under Hartley cas
                var Cg_a = Cg[a]; var Sg_a = Sg[a];
                for (int q = 0; q < L; q++)
                {
                    // xr
                    float re = Cx[q] * Cg_a[q] + Sx[q] * Sg_a[q];
                    float im = Cx[q] * Sg_a[q] - Sx[q] * Cg_a[q];
                    Hy[q] = re - im;

                    // xi
                    float re_i = Cxi[q] * Cg_a[q] + Sxi[q] * Sg_a[q];
                    float im_i = Cxi[q] * Sg_a[q] - Sxi[q] * Cg_a[q];
                    Hy2[q] = re_i - im_i;
                }
                R[a] = FHT_L.Backward(Hy);   // r_a[l]
                S[a] = FHT_L.Backward(Hy2);  // s_a[l]

                // --- xr/xi with G2 (residue a + M/2) and +1 shift along L for a ≥ M/2 ---
                var Cg2_a = Cg2[a]; var Sg2_a = Sg2[a];
                bool shift1 = (a >= (M >> 1));

                for (int q = 0; q < L; q++)
                {
                    // rotate conj(G2) by e^{-i 2π q / L} if shift1
                    float Cr = Cg2_a[q], Sr = Sg2_a[q];
                    if (shift1)
                    {
                        float cx = cL[q], s = sL[q];
                        float Cr2 = Cr * cx + Sr * s;
                        float Sr2 = -Cr * s + Sr * cx;
                        Cr = Cr2; Sr = Sr2;
                    }

                    // xr * conj(G2)
                    float re2 = Cx[q] * Cr + Sx[q] * Sr;
                    float im2 = Cx[q] * Sr - Sx[q] * Cr;
                    Hy[q] = re2 - im2;

                    // xi * conj(G2)
                    float re2i = Cxi[q] * Cr + Sxi[q] * Sr;
                    float im2i = Cxi[q] * Sr - Sxi[q] * Cr;
                    Hy2[q] = re2i - im2i;
                }
                R2[a] = FHT_L.Backward(Hy);   // r2_a[l]
                S2[a] = FHT_L.Backward(Hy2);  // s2_a[l]
            }

            // 2) For each l: FHT along a (length M) + quarter-period mixing in k
            var c = new float[2 * N];

            var vec = new float[M];
            var rev = new float[M];

            var CR = new float[M]; var SR = new float[M]; // cos/sin parts of R
            var CS = new float[M]; var SS = new float[M]; // cos/sin parts of S

            // Mixing scale along M (kept to match your matrix tests)
            float scale = 1.0f / L;

            for (int l = 0; l < L; l++)
            {
                // cos/sin for R(:, l)
                for (int a = 0; a < M; a++) vec[a] = R[a][l];
                Array.Copy(vec, rev, M); CyclicReverseInPlace(rev);
                var H1 = FHT_M.Forward(vec);
                var H1R = FHT_M.Forward(rev);
                for (int k = 0; k < M; k++) { CR[k] = 0.5f * (H1[k] + H1R[k]); SR[k] = 0.5f * (H1[k] - H1R[k]); }

                // cos/sin for S(:, l)
                for (int a = 0; a < M; a++) vec[a] = S[a][l];
                Array.Copy(vec, rev, M); CyclicReverseInPlace(rev);
                H1 = FHT_M.Forward(vec);
                H1R = FHT_M.Forward(rev);
                for (int k = 0; k < M; k++) { CS[k] = 0.5f * (H1[k] + H1R[k]); SS[k] = 0.5f * (H1[k] - H1R[k]); }

                // Branch 1: c1 = cos(φ−πk/2)·R + sin(φ−πk/2)·S
                for (int k = 0; k < M; k++)
                {
                    float re1 = c0[k] * (CR[k] + SS[k]) + s0[k] * (SR[k] - CS[k]);
                    c[l * M + k] = re1 * scale;
                }

                // cos/sin for R2(:, l) and S2(:, l)
                for (int a = 0; a < M; a++) vec[a] = R2[a][l];
                Array.Copy(vec, rev, M); CyclicReverseInPlace(rev);
                H1 = FHT_M.Forward(vec);
                H1R = FHT_M.Forward(rev);
                for (int k = 0; k < M; k++) { CR[k] = 0.5f * (H1[k] + H1R[k]); SR[k] = 0.5f * (H1[k] - H1R[k]); }

                for (int a = 0; a < M; a++) vec[a] = S2[a][l];
                Array.Copy(vec, rev, M); CyclicReverseInPlace(rev);
                H1 = FHT_M.Forward(vec);
                H1R = FHT_M.Forward(rev);
                for (int k = 0; k < M; k++) { CS[k] = 0.5f * (H1[k] + H1R[k]); SS[k] = 0.5f * (H1[k] - H1R[k]); }

                // Branch 2: c2 = −sin(φ−πk/2)·R2 + cos(φ−πk/2)·S2
                for (int k = 0; k < M; k++)
                {
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
        ///    This uses two FHTs (vector and cyclic-reverse) on weighted sequences (c0·c, s0·c).
        /// 2) For each residue a, do two L-convolutions in the Hartley domain:
        ///    xr: y1 = g_a * A1 ; y2 = g_{a+M/2}(+1 shift if a≥M/2) * (−B2)
        ///    xi: z1 = g_a * A2 ; z2 = g_{a+M/2}(+1 shift if a≥M/2) * (+B1)
        ///    Convolution in Hartley frequency uses:
        ///      H = (Cx*Cg − Sx*Sg) + (Cx*Sg + Sx*Cg)  (i.e., Re + Im in cas form)
        ///    The +1 shift along L for the G2 branch is implemented as multiplying
        ///    the window spectrum by e^{-i 2π q / L} (note the minus sign).
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
            var P = new float[M]; var PR = new float[M]; // c0-weighted
            var Q = new float[M]; var QR = new float[M]; // s0-weighted

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
                Array.Copy(P, PR, M); CyclicReverseInPlace(PR);
                var HP = FHT_M.Forward(P);
                var HPR = FHT_M.Forward(PR);
                Array.Copy(Q, QR, M); CyclicReverseInPlace(QR);
                var HQ = FHT_M.Forward(Q);
                var HQR = FHT_M.Forward(QR);
                for (int a = 0; a < M; a++)
                {
                    float cosP = 0.5f * (HP[a] + HPR[a]);
                    float sinQ = 0.5f * (HQ[a] - HQR[a]);
                    A1[a][l] = cosP + sinQ;
                }

                // --- A2 = sin_sum(c0*C1) − cos_sum(s0*C1) ---
                float[] HS = HQ, HSR = HQR; // for s0*C1
                float[] HC = HP, HCR = HPR; // for c0*C1
                for (int a = 0; a < M; a++)
                {
                    float cosS = 0.5f * (HS[a] + HSR[a]);
                    float sinC = 0.5f * (HC[a] - HCR[a]);
                    A2[a][l] = sinC - cosS;
                }

                // --- B1 = cos_sum(c0*C2) + sin_sum(s0*C2) ---
                for (int k = 0; k < M; k++) { P[k] = c0[k] * C2[k]; Q[k] = s0[k] * C2[k]; }
                Array.Copy(P, PR, M); CyclicReverseInPlace(PR);
                HP = FHT_M.Forward(P);
                HPR = FHT_M.Forward(PR);
                Array.Copy(Q, QR, M); CyclicReverseInPlace(QR);
                HQ = FHT_M.Forward(Q);
                HQR = FHT_M.Forward(QR);
                for (int a = 0; a < M; a++)
                {
                    float cosP = 0.5f * (HP[a] + HPR[a]);
                    float sinQ = 0.5f * (HQ[a] - HQR[a]);
                    B1[a][l] = cosP + sinQ;
                }

                // --- B2 = sin_sum(c0*C2) − cos_sum(s0*C2) ---
                HS = HQ; HSR = HQR;  // for s0*C2
                HC = HP; HCR = HPR;  // for c0*C2
                for (int a = 0; a < M; a++)
                {
                    float cosS = 0.5f * (HS[a] + HSR[a]);
                    float sinC = 0.5f * (HC[a] - HCR[a]);
                    B2[a][l] = sinC - cosS;
                }
            }

            // === Step 2: two L-convolutions per residue a, then deinterleave to n = a + b*M ===
            var col = new float[L];
            var colR = new float[L];
            var Cx = new float[L];
            var Sx = new float[L];
            var Hy = new float[L];

            for (int a = 0; a < M; a++)
            {
                // ----- xr: y1 = g_a * A1 -----
                Array.Copy(A1[a], col, L);
                Array.Copy(col, colR, L); CyclicReverseInPlace(colR);
                var Hx = FHT_L.Forward(col);
                var HxR = FHT_L.Forward(colR);
                for (int q = 0; q < L; q++) { Cx[q] = 0.5f * (Hx[q] + HxR[q]); Sx[q] = 0.5f * (Hx[q] - HxR[q]); }

                var Cg_a = Cg[a]; var Sg_a = Sg[a];
                for (int q = 0; q < L; q++)
                {
                    // Convolution in Hartley domain: (Cx + iSx) * (Cg + iSg)
                    float re = Cx[q] * Cg_a[q] - Sx[q] * Sg_a[q];
                    float im = Cx[q] * Sg_a[q] + Sx[q] * Cg_a[q];
                    Hy[q] = re + im; // Hartley spectrum (cas)
                }
                var y1 = FHT_L.Backward(Hy);

                // ----- xr: y2 = g_{a+M/2}(+1 shift if a≥M/2) * (−B2) -----
                Array.Copy(B2[a], col, L);
                for (int i = 0; i < L; i++) col[i] = -col[i]; // minus sign from the branch-2 formula
                Array.Copy(col, colR, L); CyclicReverseInPlace(colR);
                Hx = FHT_L.Forward(col);
                HxR = FHT_L.Forward(colR);
                for (int q = 0; q < L; q++) { Cx[q] = 0.5f * (Hx[q] + HxR[q]); Sx[q] = 0.5f * (Hx[q] - HxR[q]); }

                var Cg2_a = Cg2[a]; var Sg2_a = Sg2[a];
                bool shift1 = (a >= (M >> 1));
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
                Array.Copy(col, colR, L); CyclicReverseInPlace(colR);
                Hx = FHT_L.Forward(col);
                HxR = FHT_L.Forward(colR);
                for (int q = 0; q < L; q++) { Cx[q] = 0.5f * (Hx[q] + HxR[q]); Sx[q] = 0.5f * (Hx[q] - HxR[q]); }
                for (int q = 0; q < L; q++)
                {
                    float re = Cx[q] * Cg_a[q] - Sx[q] * Sg_a[q];
                    float im = Cx[q] * Sg_a[q] + Sx[q] * Cg_a[q];
                    Hy[q] = re + im;
                }
                var z1 = FHT_L.Backward(Hy);

                // ----- xi: z2 = g_{a+M/2}(+1 shift if a≥M/2) * (+B1) -----
                Array.Copy(B1[a], col, L);
                Array.Copy(col, colR, L); CyclicReverseInPlace(colR);
                Hx = FHT_L.Forward(col);
                HxR = FHT_L.Forward(colR);
                for (int q = 0; q < L; q++) { Cx[q] = 0.5f * (Hx[q] + HxR[q]); Sx[q] = 0.5f * (Hx[q] - HxR[q]); }
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
