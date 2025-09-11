using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Fast real 1D Weyl–Heisenberg (Gabor) transform via Hartley (FHT), matching the 2N×2N matrix (analysis exact).
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
    public class FastRealWeylHeisenbergTransformHartley : TransformBaseFloat, ITransform
    {
        #region Private data
        private readonly FastHartleyTransform FHT_L; // along L
        private readonly FastHartleyTransform FHT_M; // along M

        private float[] g0;   // prototype window, length N
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

        private Direction direction = Direction.Vertical;
        #endregion

        #region Ctor / Properties
        public FastRealWeylHeisenbergTransformHartley(float[] g0, int M, bool normalized = false, Direction direction = Direction.Vertical)
        {
            if (g0 == null) throw new ArgumentNullException(nameof(g0));
            if (M <= 0 || (M & 1) != 0) throw new ArgumentException("M must be positive and even.", nameof(M));

            this.N = g0.Length;
            if (N % M != 0) throw new ArgumentException("g0.Length must be divisible by M (N = M * L).", nameof(g0));

            this.M = M;
            this.L = N / M;
            this.g0 = (float[])g0.Clone();
            this.direction = direction;

            // Hartley transforms
            this.FHT_L = new FastHartleyTransform(normalized, direction);
            this.FHT_M = new FastHartleyTransform(normalized, direction);

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

        public bool Normalized
        {
            get => FHT_L.Normalized && FHT_M.Normalized;
            set { FHT_L.Normalized = value; FHT_M.Normalized = value; }
        }

        public override Direction Direction
        {
            get => direction;
            set => direction = value;
        }

        public float[] Window => (float[])g0.Clone();

        public void SetWindow(float[] newWindow)
        {
            if (newWindow == null) throw new ArgumentNullException(nameof(newWindow));
            if (newWindow.Length != N) throw new ArgumentException("Window length mismatch.", nameof(newWindow));
            g0 = (float[])newWindow.Clone();
            PrecomputeWindowSpectra();
        }
        #endregion

        #region Helpers
        private static void CyclicReverseInPlace(float[] v)
        {
            int n = v.Length, half = n / 2;
            for (int i = 1; i <= half; i++)
            {
                int a = i, b = n - i;
                var t = v[a]; v[a] = v[b]; v[b] = t;
            }
        }

        private void PrecomputeWindowSpectra()
        {
            var Ga = new float[L];
            var G2a = new float[L];
            var tmp = new float[L];

            for (int a = 0; a < M; a++)
            {
                // residue a
                for (int b = 0; b < L; b++) Ga[b] = g0[a + b * M];

                var H = FHT_L.Forward(Ga);
                Array.Copy(Ga, tmp, L); CyclicReverseInPlace(tmp);
                var HR = FHT_L.Forward(tmp);

                var C = new float[L];
                var S = new float[L];

                for (int q = 0; q < L; q++)
                {
                    C[q] = 0.5f * (H[q] + HR[q]); // Re
                    S[q] = 0.5f * (H[q] - HR[q]); // Im-like (sign matched for conj)
                }

                Cg[a] = C; Sg[a] = S;

                // residue a' = a + M/2
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

        #region Forward (analysis) — y[2N] -> c[2N]
        /// <summary>
        /// Forward (analysis): y[2N] -> c[2N] equals G^T * y for your 2N×2N matrix.
        /// If y.Length == N, the second half is assumed zero (xi = 0).
        /// </summary>
        public override float[] Forward(float[] y)
        {
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (y.Length != N && y.Length != 2 * N)
                throw new ArgumentException($"Input length must be N={N} or 2N={2 * N}.");

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
            var Xrv = new float[L];
            var Hy = new float[L];
            var Hy2 = new float[L];

            for (int a = 0; a < M; a++)
            {
                // --- xr with Ga ---
                for (int b = 0; b < L; b++) Xa[b] = xr[a + b * M];

                var Hx = FHT_L.Forward(Xa);
                Array.Copy(Xa, Xrv, L); CyclicReverseInPlace(Xrv);
                var HxR = FHT_L.Forward(Xrv);

                var Cx = new float[L];
                var Sx = new float[L];

                for (int q = 0; q < L; q++)
                {
                    Cx[q] = 0.5f * (Hx[q] + HxR[q]);
                    Sx[q] = 0.5f * (Hx[q] - HxR[q]);
                }

                // --- xi with Ga ---
                for (int b = 0; b < L; b++) Xa[b] = xi[a + b * M];

                var Hi = FHT_L.Forward(Xa);
                Array.Copy(Xa, Xrv, L); CyclicReverseInPlace(Xrv);
                var HiR = FHT_L.Forward(Xrv);

                var Cxi = new float[L];
                var Sxi = new float[L];

                for (int q = 0; q < L; q++)
                {
                    Cxi[q] = 0.5f * (Hi[q] + HiR[q]);
                    Sxi[q] = 0.5f * (Hi[q] - HiR[q]);
                }

                // Y = X * conj(G) for Ga
                var Cg_a = Cg[a]; var Sg_a = Sg[a];

                for (int q = 0; q < L; q++)
                {
                    float re = Cx[q] * Cg_a[q] + Sx[q] * Sg_a[q];
                    float im = Cx[q] * Sg_a[q] - Sx[q] * Cg_a[q];
                    Hy[q] = re - im; // Hartley spectrum for correlation of xr with Ga

                    float re_i = Cxi[q] * Cg_a[q] + Sxi[q] * Sg_a[q];
                    float im_i = Cxi[q] * Sg_a[q] - Sxi[q] * Cg_a[q];
                    Hy2[q] = re_i - im_i; // Hartley spectrum for correlation of xi with Ga
                }

                R[a] = FHT_L.Backward(Hy);     // r_a[l]
                S[a] = FHT_L.Backward(Hy2);    // s_a[l]

                // Y = X * conj(G) for Ga2 (residue a+M/2), with extra +1 shift along L when a ≥ M/2
                var Cg2_a = Cg2[a]; var Sg2_a = Sg2[a];
                bool shift1 = (a >= (M >> 1));

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

                R2[a] = FHT_L.Backward(Hy);    // r2_a[l]
                S2[a] = FHT_L.Backward(Hy2);   // s2_a[l]
            }

            // 2) For each l: FHT along a (length M) + quarter-period mixing
            var c = new float[2 * N];

            var vec = new float[M];
            var rev = new float[M];

            var CR = new float[M]; var SR = new float[M];
            var CS = new float[M]; var SS = new float[M];

            // scale for M-stage mixing:
            float scale = 1.0f / L;

            for (int l = 0; l < L; l++)
            {
                // --- branch 1 (G1): from R,S ---
                for (int a = 0; a < M; a++) vec[a] = R[a][l];
                Array.Copy(vec, rev, M); CyclicReverseInPlace(rev);
                var H1 = FHT_M.Forward(vec);
                var H1R = FHT_M.Forward(rev);
                for (int k = 0; k < M; k++) { CR[k] = 0.5f * (H1[k] + H1R[k]); SR[k] = 0.5f * (H1[k] - H1R[k]); }

                for (int a = 0; a < M; a++) vec[a] = S[a][l];
                Array.Copy(vec, rev, M); CyclicReverseInPlace(rev);
                H1 = FHT_M.Forward(vec);
                H1R = FHT_M.Forward(rev);
                for (int k = 0; k < M; k++) { CS[k] = 0.5f * (H1[k] + H1R[k]); SS[k] = 0.5f * (H1[k] - H1R[k]); }

                for (int k = 0; k < M; k++)
                {
                    // c1 = cos(φ−πk/2)·R + sin(φ−πk/2)·S
                    float re1 = c0[k] * (CR[k] + SS[k]) + s0[k] * (SR[k] - CS[k]);
                    c[l * M + k] = re1 * scale;
                }

                // --- branch 2 (G2): from R2,S2) ---
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

                for (int k = 0; k < M; k++)
                {
                    // c2 = −sin(φ−πk/2)·R2 + cos(φ−πk/2)·S2
                    float re2 = c0[k] * (-SR[k] + CS[k]) + s0[k] * (CR[k] + SS[k]);
                    c[l * M + k + N] = re2 * scale;
                }
            }

            return c;
        }
        #endregion

        #region Backward (synthesis) — exact O(N^2) matrix form
        /// <summary>
        /// Backward (synthesis): c[2N] -> y[2N], exact matrix-form (O(N^2)) using Special.Cas (Hartley-only).
        /// </summary>
        public override float[] Backward(float[] c)
        {
            if (c == null) throw new ArgumentNullException(nameof(c));
            if (c.Length != 2 * N) throw new ArgumentException($"Input length must be 2N = {2 * N}.");

            var xr = new float[N];
            var xi = new float[N];

            float twoPiOverM = 2.0f * Maths.Pi / M;
            float halfPi = 0.5f * Maths.Pi;

            for (int n = 0; n < N; n++)
            {
                float sumR = 0f, sumI = 0f;

                for (int l = 0; l < L; l++)
                {
                    for (int k = 0; k < M; k++)
                    {
                        int u = l * M + k;
                        int i = Maths.Mod(n - l * M, N);
                        int j = Maths.Mod(n + (M >> 1) - l * M, N);

                        float theta = twoPiOverM * k * (n - 0.25f * M);

                        // cosθ, sinθ via cas
                        float cas1 = Special.Cas(theta);            // cosθ + sinθ
                        float cas2 = Special.Cas(theta - halfPi);   // sinθ - cosθ
                        float cosT = 0.5f * (cas1 - cas2);
                        float sinT = 0.5f * (cas1 + cas2);

                        float gi = g0[i];
                        float gj = g0[j];

                        // G1 contributions
                        sumR += gi * cosT * c[u];
                        sumI += gi * sinT * c[u];

                        // G2 contributions (i * gj * e^{iθ})
                        sumR += (-gj * sinT) * c[u + N];
                        sumI += (gj * cosT) * c[u + N];
                    }
                }

                xr[n] = sumR;
                xi[n] = sumI;
            }

            var y = new float[2 * N];
            Buffer.BlockCopy(xr, 0, y, 0, sizeof(float) * N);
            Buffer.BlockCopy(xi, 0, y, sizeof(float) * N, sizeof(float) * N);
            return y;
        }
        #endregion
    }
}
