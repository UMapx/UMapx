using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines fast real Weyl-Heisenberg transform.
    /// </summary>
    /// <remarks>
    /// The class represents a computationally efficient implementation of one-dimensional and two-dimensional discrete real orthogonal
    /// Weyl-Heisenberg transforms. This implementation was designed and developed by Valery Asiryan, Yerevan, Armenia (2025).
    /// </remarks>
    public partial class FastRealWeylHeisenbergTransform
    {
        #region Internal methods and classes

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
        /// 2) Along M (frequency shifts): FHT + reverse trick to get cos/sin for each of r,s,
        ///    then quarter-period mixing using k mod 4 via c0[k]=cos(πk/2), s0[k]=sin(πk/2) ∈ {0,±1}.
        /// 3) Branch formulas:
        ///      c1 =  cos(φ − πk/2)·R  +  sin(φ − πk/2)·S
        ///      c2 = −sin(φ − πk/2)·R2 +  cos(φ − πk/2)·S2
        /// 4) A final scale = 1/L is applied to match your matrix convention.
        /// </remarks>
        public float[] FRWHT(float[] y, RealPolyphaseCache cache)
        {
            var N = cache.N;
            var L = cache.L;
            var M = cache.M;

            // quarter-period tables
            var c0 = new sbyte[M];
            var s0 = new sbyte[M];

            for (int k = 0; k < M; k++)
            {
                int r = k & 3;
                c0[k] = (sbyte)((r == 0) ? +1 : (r == 2) ? -1 : 0);
                s0[k] = (sbyte)((r == 1) ? +1 : (r == 3) ? -1 : 0);
            }

            // precompute cos/sin for e^{-i 2π q/L}
            var cL = new float[L];
            var sL = new float[L];

            for (int q = 0; q < L; q++)
            {
                double ang = 2.0 * Math.PI * q / L;
                cL[q] = (float)Math.Cos(ang);
                sL[q] = (float)Math.Sin(ang);
            }

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
                var Hx = FHT.Forward(Xa);

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

                var Hi = FHT.Forward(Xa);

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
                var Cg_a = cache.Cg[a]; var Sg_a = cache.Sg[a];

                for (int q = 0; q < L; q++)
                {
                    float re = Cx[q] * Cg_a[q] + Sx[q] * Sg_a[q];
                    float im = Cx[q] * Sg_a[q] - Sx[q] * Cg_a[q];
                    Hy[q] = re - im;   // xr with Ga

                    float re_i = Cxi[q] * Cg_a[q] + Sxi[q] * Sg_a[q];
                    float im_i = Cxi[q] * Sg_a[q] - Sxi[q] * Cg_a[q];
                    Hy2[q] = re_i - im_i; // xi with Ga
                }

                R[a] = FHT.Backward(Hy);   // r_a[l]
                S[a] = FHT.Backward(Hy2);  // s_a[l]

                // Y = X * conj(G2) for Ga2 (residue a+M/2), with extra +1 shift along L when a ≥ M/2
                var Cg2_a = cache.Cg2[a]; var Sg2_a = cache.Sg2[a];
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

                R2[a] = FHT.Backward(Hy);   // r2_a[l]
                S2[a] = FHT.Backward(Hy2);  // s2_a[l]
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
                var H1 = FHT.Forward(vec);
                for (int k = 0; k < M; k++)
                {
                    int krev = (k == 0) ? 0 : (M - k);
                    float HR = H1[krev];
                    CR[k] = 0.5f * (H1[k] + HR);
                    SR[k] = 0.5f * (H1[k] - HR);
                }

                for (int a = 0; a < M; a++) vec[a] = S[a][l];
                H1 = FHT.Forward(vec);
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
                H1 = FHT.Forward(vec);
                for (int k = 0; k < M; k++)
                {
                    int krev = (k == 0) ? 0 : (M - k);
                    float HR = H1[krev];
                    CR[k] = 0.5f * (H1[k] + HR);
                    SR[k] = 0.5f * (H1[k] - HR);
                }

                for (int a = 0; a < M; a++) vec[a] = S2[a][l];
                H1 = FHT.Forward(vec);
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
        /// 1) For each l, recover four k-sums via 2×FHT with quarter-period tables:
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
        public float[] IFRWHT(float[] c, RealPolyphaseCache cache)
        {
            if (c == null) throw new ArgumentNullException(nameof(c));

            var N = cache.N;
            var L = cache.L;
            var M = cache.M;

            // quarter-period tables
            var c0 = new sbyte[M];
            var s0 = new sbyte[M];

            for (int k = 0; k < M; k++)
            {
                int r = k & 3;
                c0[k] = (sbyte)((r == 0) ? +1 : (r == 2) ? -1 : 0);
                s0[k] = (sbyte)((r == 1) ? +1 : (r == 3) ? -1 : 0);
            }

            // precompute cos/sin for e^{-i 2π q/L}
            var cL = new float[L];
            var sL = new float[L];

            for (int q = 0; q < L; q++)
            {
                double ang = 2.0 * Math.PI * q / L;
                cL[q] = (float)Math.Cos(ang);
                sL[q] = (float)Math.Sin(ang);
            }

            var xr = new float[N];
            var xi = new float[N];

            // Undo the Forward mixing scale (Forward stored re*scale, here we multiply by 1/scale = L).
            float invScaleM = 1.0f / L;

            // === Step 1: recover A1, A2, B1, B2 for each l via two FHT calls ===
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

                var HP = FHT.Forward(P);
                var HQ = FHT.Forward(Q);

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

                HP = FHT.Forward(P);
                HQ = FHT.Forward(Q);

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
                var Hx = FHT.Forward(col);
                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);
                    float HR = Hx[qrev];
                    Cx[q] = 0.5f * (Hx[q] + HR);
                    Sx[q] = 0.5f * (Hx[q] - HR);
                }

                var Cg_a = cache.Cg[a]; var Sg_a = cache.Sg[a];
                for (int q = 0; q < L; q++)
                {
                    float re = Cx[q] * Cg_a[q] - Sx[q] * Sg_a[q];
                    float im = Cx[q] * Sg_a[q] + Sx[q] * Cg_a[q];
                    Hy[q] = re + im; // Hartley spectrum (cas)
                }
                var y1 = FHT.Backward(Hy);

                // ----- xr: y2 = g_{a+M/2}(+1 shift if a≥M/2) * (−B2) -----
                Array.Copy(B2[a], col, L);
                for (int i = 0; i < L; i++) col[i] = -col[i]; // minus sign from the branch-2 formula
                Hx = FHT.Forward(col);
                for (int q = 0; q < L; q++)
                {
                    int qrev = (q == 0) ? 0 : (L - q);
                    float HR = Hx[qrev];
                    Cx[q] = 0.5f * (Hx[q] + HR);
                    Sx[q] = 0.5f * (Hx[q] - HR);
                }

                var Cg2_a = cache.Cg2[a]; var Sg2_a = cache.Sg2[a];
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
                var y2 = FHT.Backward(Hy);

                // Scatter back to xr at n = a + b*M
                for (int b = 0; b < L; b++)
                {
                    int n = a + b * M;
                    xr[n] = y1[b] + y2[b];
                }

                // ----- xi: z1 = g_a * A2 -----
                Array.Copy(A2[a], col, L);
                Hx = FHT.Forward(col);
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
                var z1 = FHT.Backward(Hy);

                // ----- xi: z2 = g_{a+M/2}(+1 shift if a≥M/2) * (+B1) -----
                Array.Copy(B1[a], col, L);
                Hx = FHT.Forward(col);
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
                var z2 = FHT.Backward(Hy);

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

        /// <summary>
        /// Polyphase/Hartley cache for the fast real Weyl–Heisenberg (Gabor) transform.
        /// The cache stores, for every residue a, the Hartley-domain
        /// spectra of the window's polyphase components along the L-dimension.
        /// 
        /// Notation:
        ///   N = total length, M = number of frequency shifts (must be even), L = N / M
        ///   n = a + b·M  (polyphase indexing: residue a and time index b)
        /// 
        /// What is cached (per residue a, per frequency bin q in 0..L-1):
        ///   G_a(q)  = Cg[a][q]  + i · Sg[a][q]
        ///   G2_a(q) = Cg2[a][q] + i · Sg2[a][q], where a' = (a + M/2) mod M
        /// 
        /// The pair (C,S) is a "complex-like" decomposition in the Hartley (cas) domain:
        /// we recover cosine/sine sums via one FHT of the polyphase column and its
        /// cyclic-reversal using the identity H_rev[q] = H[(L - q) mod L]:
        ///   C[q] = 0.5 · (H[q] + H_rev[q])  == Σ x[b] · cos(2π·q·b/L)
        ///   S[q] = 0.5 · (H[q] - H_rev[q])  == Σ x[b] · sin(2π·q·b/L)
        /// 
        /// These cached spectra are used in both Forward and Backward to avoid repeated
        /// window transforms and to keep numerics stable/deterministic.
        /// </summary>
        public sealed class RealPolyphaseCache
        {
            /// <summary>
            /// Total signal length N. Must satisfy N = M x L.
            /// </summary>
            public readonly int N;

            /// <summary>
            /// Number of frequency shifts (must be even).
            /// </summary>
            public readonly int M;

            /// <summary>
            /// Number of time shifts (L = N / M).
            /// </summary>
            public readonly int L;

            // -------------------- Cached spectra --------------------

            /// <summary>
            /// Real (cosine) part of the Hartley "spectrum" G_a(q) for residue a (main branch).
            /// Size: M arrays of length L.
            /// </summary>
            public readonly float[][] Cg;

            /// <summary>
            /// Imag-like (sine) part of the Hartley "spectrum" G_a(q) for residue a (main branch).
            /// The sign convention is chosen so that multiplication by conj(G) in the
            /// Hartley domain maps to the (re − im) combination used in Forward().
            /// Size: M arrays of length L.
            /// </summary>
            public readonly float[][] Sg;

            /// <summary>
            /// Real (cosine) part of the Hartley "spectrum" G2_a(q) for residue a (half-shifted branch).
            /// Here G2 corresponds to residue a' = a + M/2 (mod M).
            /// Size: M arrays of length L.
            /// </summary>
            public readonly float[][] Cg2;

            /// <summary>
            /// Imag-like (sine) part of the Hartley "spectrum" G2_a(q) for residue a (half-shifted branch).
            /// Size: M arrays of length L.
            /// </summary>
            public readonly float[][] Sg2;

            /// <summary>
            /// Construct a cache instance from precomputed arrays. Typically produced by <see cref="Build"/>.
            /// </summary>
            /// <param name="N">Total signal length (N = M·L).</param>
            /// <param name="M">Number of frequency shifts (even).</param>
            /// <param name="L">Number of time shifts (L = N / M).</param>
            /// <param name="Cg">Real parts of G_a(q).</param>
            /// <param name="Sg">Imag-like parts of G_a(q).</param>
            /// <param name="Cg2">Real parts of G2_a(q) for a' = a + M/2.</param>
            /// <param name="Sg2">Imag-like parts of G2_a(q) for a' = a + M/2.</param>
            public RealPolyphaseCache(int N, int M, int L, float[][] Cg, float[][] Sg, float[][] Cg2, float[][] Sg2)
            {
                this.N = N;
                this.M = M;
                this.L = L;
                this.Cg = Cg;
                this.Sg = Sg;
                this.Cg2 = Cg2;
                this.Sg2 = Sg2;
            }

            /// <summary>
            /// Build a <see cref="RealPolyphaseCache"/> from a window function and grid (N, M).
            /// 
            /// Steps:
            ///  1) Generate the length-N analysis window g0 using the same factory
            ///     as the complex implementation.
            ///  2) Orthonormalize g0 via Zak-domain orthogonalization to obtain WH-orthonormal g.
            ///     (This ensures the real/Hartley path uses an identically normalized window.)
            ///  3) For each residue a in 0..M-1:
            ///       • Form the polyphase column Ga[b] = g[a + b·M], b = 0..L-1
            ///       • Compute its length-L FHT: H = FHT(Ga)
            ///       • Recover cosine/sine sums via mirror:
            ///           C[q] = 0.5·(H[q] + H[(L − q) mod L])
            ///           S[q] = 0.5·(H[q] − H[(L − q) mod L])
            ///       • Store as Cg[a], Sg[a]
            ///     Then repeat for the half-shift residue a' = a + M/2 (mod M) to get Cg2[a], Sg2[a].
            /// 
            /// Notes:
            ///  • M must be even so that a' is well-defined.
            ///  • We deliberately use the orthonormalized window g for both G and G2 branches.
            ///    Mixing g and g0 here would break orthogonality and lead to amplitude mismatches.
            ///  • This method uses the outer transform's single FHT instance (shared),
            ///    which is safe because nested types in C# can access private members of the enclosing type.
            /// </summary>
            /// <param name="N">Total signal length (must be divisible by Mloc).</param>
            /// <param name="Mloc">Number of frequency shifts M (must be even).</param>
            /// <param name="window">Window function (e.g., Gaussian, Hann, etc.).</param>
            /// <returns>Initialized polyphase/Hartley cache.</returns>
            /// <exception cref="ArgumentException">Thrown if N is not divisible by Mloc or Mloc is not even.</exception>
            public static RealPolyphaseCache Build(int N, int Mloc, IWindow window)
            {
                // Validate grid: N must factor as Mloc·L and Mloc must be even.
                int L = N / Mloc;
                if (L * Mloc != N) throw new ArgumentException("N must be divisible by M");
                if ((Mloc & 1) != 0) throw new ArgumentException("M must be even");

                // 1) Build the prototype analysis window g0 (length N) using the same
                //    factory as the complex WH transform to keep normalization consistent.
                var g0 = FastWeylHeisenbergTransform.Packet(window, N);

                // 2) Zak-domain orthogonalization: produce WH-orthonormal window g.
                //    This matches the "Matrix(..., true)" behavior in the slow reference path.
                var zakOrth = new FastZakTransform(Mloc);
                var g = zakOrth.Orthogonalize(g0); // real-valued, length N

                // Temporary buffers for polyphase columns (length L).
                var Ga = new float[L];
                var G2a = new float[L];

                // Allocate output arrays.
                var Cg = new float[Mloc][];
                var Sg = new float[Mloc][];
                var Cg2 = new float[Mloc][];
                var Sg2 = new float[Mloc][];

                for (int a = 0; a < Mloc; a++)
                {
                    // -------- Main branch: residue a --------
                    // Build polyphase column Ga[b] = g[a + b·M], b = 0..L-1
                    for (int b = 0; b < L; b++) Ga[b] = g[a + b * Mloc];

                    // One length-L FHT; mirrored bins give us cos/sin projections.
                    var H = FHT.Forward(Ga);

                    var C = new float[L];
                    var S = new float[L];

                    // Recover cos/sin sums from Hartley bins and their "reverse":
                    //   H_rev[q] = H[(L − q) mod L]
                    //   C[q] = 0.5·(H[q] + H_rev[q])  → Σ Ga[b]·cos(...)
                    //   S[q] = 0.5·(H[q] − H_rev[q])  → Σ Ga[b]·sin(...)
                    for (int q = 0; q < L; q++)
                    {
                        int qrev = (q == 0) ? 0 : (L - q); // maps 0→0, (L/2)→(L/2) when L is even
                        float HR = H[qrev];
                        C[q] = 0.5f * (H[q] + HR);
                        S[q] = 0.5f * (H[q] - HR);
                    }

                    Cg[a] = C;
                    Sg[a] = S;

                    // -------- Half-shift branch: residue a' = a + M/2 (mod M) --------
                    int ap = (a + (Mloc >> 1)) % Mloc;

                    // IMPORTANT: use the orthonormalized window g here as well
                    // (not g0), otherwise branches G and G2 would be inconsistent.
                    for (int b = 0; b < L; b++) G2a[b] = g[ap + b * Mloc];

                    var H2 = FHT.Forward(G2a);

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

                // Return the completed cache object.
                return new RealPolyphaseCache(N, Mloc, L, Cg, Sg, Cg2, Sg2);
            }
        }

        /// <summary>
        /// UMapx fast Hartley transform.
        /// </summary>
        private static readonly FastHartleyTransform FHT = new FastHartleyTransform(false, Direction.Vertical);

        #endregion
    }
}
