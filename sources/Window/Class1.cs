using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Fast real 1D Weyl–Heisenberg (Gabor) transform via Hartley (FHT), matching the 2N×2N matrix.
    /// </summary>
    /// <remarks>
    /// Factorization (analysis): polyphase split N = M * L; correlation along L via Hartley-DFT (two FHT_L calls);
    /// then FHT_M along a and quarter-period (-π/2) mixing via {cos(πk/2), sin(πk/2)} ∈ {0,±1}.
    /// Second branch (M/2 shift) is handled via residues (a + M/2) and an extra (-1)^k.
    /// Synthesis here is an exact O(N^2) reference (matrix-accurate) using Special.Cas; no complex arithmetic.
    /// Uses only your FastHartleyTransform for all DFT-like steps.
    /// </remarks>
    [Serializable]
    public class FastRealWeylHeisenbergTransformHartley : TransformBaseFloat, ITransform
    {
        #region Private data
        private readonly FastHartleyTransform FHT_L;   // Hartley along L
        private readonly FastHartleyTransform FHT_M;   // Hartley along M
        private float[] g0;                            // prototype window, length N
        private readonly int M;                        // number of frequency shifts (even)
        private readonly int L;                        // time steps, N = M * L
        private readonly int N;                        // signal length

        // quarter-period twiddles per k: cos(πk/2), sin(πk/2) ∈ {0,±1}
        private readonly sbyte[] c0;
        private readonly sbyte[] s0;

        // half-period twiddle (-1)^k for M/2 branch
        private readonly float[] s2;

        // conj DFT of window polyphases along L, for each residue a
        private readonly float[][] Cg;    // Re spectrum for G_a
        private readonly float[][] Sg;    // Im-like component for G_a (sign adjusted for Hartley usage)
        private readonly float[][] Cg2;   // for residue a' = a + M/2
        private readonly float[][] Sg2;

        private readonly float eps;       // guard (used in synthesis fallbacks only)
        private Direction direction = Direction.Vertical;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast real 1D WH transform (2N↔2N) using Hartley transforms only.
        /// </summary>
        /// <param name="g0">Prototype window of length N (N must be divisible by M)</param>
        /// <param name="M">Number of frequency shifts (must be even)</param>
        /// <param name="normalized">Use orthonormal Hartley scaling (true = 1/√n forward/backward)</param>
        /// <param name="eps">Small guard for divisions in synthesis (default 1e-7)</param>
        /// <param name="direction">Processing direction flag (kept for API symmetry)</param>
        public FastRealWeylHeisenbergTransformHartley(float[] g0, int M, bool normalized = true, float eps = 1e-7f, Direction direction = Direction.Vertical)
        {
            if (g0 == null) throw new ArgumentNullException(nameof(g0));
            if (M <= 0) throw new ArgumentException("M must be positive.", nameof(M));
            if ((M & 1) != 0) throw new ArgumentException("M must be even (M % 2 == 0).", nameof(M));

            this.N = g0.Length;
            if (N % M != 0) throw new ArgumentException("g0.Length must be divisible by M (N = M * L).", nameof(g0));

            this.M = M;
            this.L = N / M;
            this.g0 = (float[])g0.Clone();
            this.eps = Math.Max(0.0f, eps);
            this.direction = direction;

            // Hartley transforms (your implementation)
            this.FHT_L = new FastHartleyTransform(normalized, direction);
            this.FHT_M = new FastHartleyTransform(normalized, direction);

            // quarter- and half-period twiddles
            this.c0 = new sbyte[M];
            this.s0 = new sbyte[M];
            this.s2 = new float[M];

            for (int k = 0; k < M; k++)
            {
                int r = k & 3; // k mod 4
                c0[k] = (sbyte)((r == 0) ? +1 : (r == 2) ? -1 : 0);
                s0[k] = (sbyte)((r == 1) ? +1 : (r == 3) ? -1 : 0);
                s2[k] = ((k & 1) == 0) ? +1f : -1f; // (-1)^k
            }

            // allocate spectra
            this.Cg = new float[M][];
            this.Sg = new float[M][];
            this.Cg2 = new float[M][];
            this.Sg2 = new float[M][];

            // precompute window spectra along L for each residue a (and a+M/2)
            PrecomputeWindowSpectra();
        }

        /// <summary>Whether the internal Hartley transforms are normalized.</summary>
        public bool Normalized
        {
            get => FHT_L.Normalized && FHT_M.Normalized;
            set { FHT_L.Normalized = value; FHT_M.Normalized = value; }
        }

        /// <summary>Processing direction (kept for API symmetry).</summary>
        public override Direction Direction
        {
            get => direction;
            set => direction = value;
        }

        /// <summary>Returns a copy of the current window.</summary>
        public float[] Window => (float[])g0.Clone();

        /// <summary>Updates the window and refreshes all precomputed spectra.</summary>
        public void SetWindow(float[] newWindow)
        {
            if (newWindow == null) throw new ArgumentNullException(nameof(newWindow));
            if (newWindow.Length != N) throw new ArgumentException("New window must have length N.", nameof(newWindow));
            this.g0 = (float[])newWindow.Clone();
            PrecomputeWindowSpectra();
        }

        private static void CyclicReverseInPlace(float[] v)
        {
            int n = v.Length;
            int half = n / 2;
            for (int i = 1; i <= half; i++)
            {
                int a = i, b = n - i;
                float t = v[a]; v[a] = v[b]; v[b] = t;
            }
        }

        /// <summary>
        /// Precomputes conj DFT along L for each residue a: Ga[b] = g0[a + b*M], b=0..L-1.
        /// Two-Hartley trick to recover (C,S); also for the M/2-shifted residue.
        /// </summary>
        private void PrecomputeWindowSpectra()
        {
            var Ga = new float[L];
            var G2a = new float[L];
            var tmp = new float[L];

            for (int a = 0; a < M; a++)
            {
                // --- Residue a ---
                for (int b = 0; b < L; b++) Ga[b] = g0[a + b * M];

                var H = FHT_L.Forward(Ga);         // H = C + S
                Array.Copy(Ga, tmp, L);
                CyclicReverseInPlace(tmp);
                var HR = FHT_L.Forward(tmp);        // HR = C - S

                var C = new float[L];
                var S = new float[L];
                for (int q = 0; q < L; q++)
                {
                    C[q] = 0.5f * (H[q] + HR[q]);  // Re
                    S[q] = 0.5f * (H[q] - HR[q]);  // Im-like (sign matched for conj usage)
                }
                Cg[a] = C; Sg[a] = S;

                // --- Residue a' = a + M/2 (wrapped) ---
                int ap = (a + (M >> 1)) % M;
                for (int b = 0; b < L; b++) G2a[b] = g0[ap + b * M];

                var H2 = FHT_L.Forward(G2a);
                Array.Copy(G2a, tmp, L);
                CyclicReverseInPlace(tmp);
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

        #region Forward (analysis) 2N -> 2N
        /// <summary>
        /// Forward (analysis): y[2N] -> c[2N] equals G^T * y for your 2N×2N matrix.
        /// If y.Length == N, the second half is assumed zero (xi = 0).
        /// </summary>
        public override float[] Forward(float[] y)
        {
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (y.Length != N && y.Length != 2 * N)
                throw new ArgumentException($"Input length must be N={N} or 2N={2 * N}.");

            // Split input into top/bottom halves (xr, xi)
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
            }

            // 1) Correlation along L for each residue a (for both branches and both halves)
            var R = new float[M][]; // from xr with Ga
            var S = new float[M][]; // from xi with Ga
            var R2 = new float[M][]; // from xr with G_{a+M/2}
            var S2 = new float[M][]; // from xi with G_{a+M/2}

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

                // Y = X * conj(G) for Ga2 (residue a+M/2)
                var Cg2_a = Cg2[a]; var Sg2_a = Sg2[a];
                for (int q = 0; q < L; q++)
                {
                    float re2 = Cx[q] * Cg2_a[q] + Sx[q] * Sg2_a[q];
                    float im2 = Cx[q] * Sg2_a[q] - Sx[q] * Cg2_a[q];
                    Hy[q] = re2 - im2;

                    float re2i = Cxi[q] * Cg2_a[q] + Sxi[q] * Sg2_a[q];
                    float im2i = Cxi[q] * Sg2_a[q] - Sxi[q] * Cg2_a[q];
                    Hy2[q] = re2i - im2i;
                }
                R2[a] = FHT_L.Backward(Hy);    // r2_a[l]
                S2[a] = FHT_L.Backward(Hy2);   // s2_a[l]
            }

            // 2) For each l: FHT along a (length M) + quarter-period mixing
            var c = new float[2 * N];

            var vec = new float[M];
            var rev = new float[M];

            var H1 = new float[M];
            var H1R = new float[M];

            var CR = new float[M]; var SR = new float[M];
            var CS = new float[M]; var SS = new float[M];

            for (int l = 0; l < L; l++)
            {
                // --- branch 1 (G1): from R,S ---
                for (int a = 0; a < M; a++) vec[a] = R[a][l];
                Array.Copy(vec, rev, M); CyclicReverseInPlace(rev);
                H1 = FHT_M.Forward(vec);
                H1R = FHT_M.Forward(rev);
                for (int k = 0; k < M; k++) { CR[k] = 0.5f * (H1[k] + H1R[k]); SR[k] = 0.5f * (H1[k] - H1R[k]); }

                for (int a = 0; a < M; a++) vec[a] = S[a][l];
                Array.Copy(vec, rev, M); CyclicReverseInPlace(rev);
                H1 = FHT_M.Forward(vec);
                H1R = FHT_M.Forward(rev);
                for (int k = 0; k < M; k++) { CS[k] = 0.5f * (H1[k] + H1R[k]); SS[k] = 0.5f * (H1[k] - H1R[k]); }

                for (int k = 0; k < M; k++)
                {
                    // c[u] = (cos⋅R + sin⋅S) with quarter shift: θ -> θ - π/2
                    // cos(θ-π/2) = c0*C + s0*S ;  sin(θ-π/2) = c0*S - s0*C
                    float re1 = c0[k] * (CR[k] + SS[k]) + s0[k] * (SR[k] - CS[k]);
                    int u = l * M + k;
                    c[u] = re1 / 2;

                }

                // --- branch 2 (G2): from R2,S2 (extra (-1)^k factor) ---
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

                // --- branch 2 (G2) ---
                for (int k = 0; k < M; k++)
                {
                    float re2 = c0[k] * (-SR[k] + CS[k]) + s0[k] * (CR[k] + SS[k]);
                    int u = l * M + k;
                    c[u + N] = re2 / 2;
                }
            }

            return c;
        }
        #endregion

        #region Backward (synthesis) 2N -> 2N (exact via matrix)
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

                        // θ(n,k) = 2π (k/M) (n - M/4)
                        float theta = twoPiOverM * k * (n - 0.25f * M);

                        // cosθ, sinθ via cas to stay Hartley-only
                        float cas1 = Special.Cas(theta);            // cosθ + sinθ
                        float cas2 = Special.Cas(theta - halfPi);   // sinθ - cosθ
                        float cosT = 0.5f * (cas1 - cas2);
                        float sinT = 0.5f * (cas1 + cas2);

                        float gi = g0[i];
                        float gj = g0[j];

                        // G1 branch contributions
                        sumR += gi * cosT * c[u];
                        sumI += gi * sinT * c[u];

                        // G2 branch contributions (i * g0[j] * e^{iθ})
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
