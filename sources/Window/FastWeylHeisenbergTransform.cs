using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines fast Weyl-Heisenberg transform.
    /// <remarks>
    /// The class represents a computationally efficient implementation of one-dimensional and two-dimensional discrete orthogonal
    /// Weyl-Heisenberg transforms. This implementation was designed and developed by Valery Asiryan, Yerevan, Armenia (2025).
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastWeylHeisenbergTransform : WeylHeisenbergTransform, IWindowTransform, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes fast Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N]</param>
        /// <param name="direction">Processing direction</param>
        public FastWeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical) : base(window, m, direction) { }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            var cache = PolyphaseCache.Build(N, this.m, this.window);
            return FWHT(A, cache);
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length / 2;
            var cache = PolyphaseCache.Build(N, this.m, this.window);
            return IFWHT(B, cache);
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public override Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0); // rows
            int M = A.GetLength(1); // cols

            PolyphaseCache cacheCols = null; // U for length N (vertical)
            PolyphaseCache cacheRows = null; // V for length M (horizontal)

            if (direction == Direction.Both || direction == Direction.Vertical)
                cacheCols = PolyphaseCache.Build(N, this.m, this.window);

            if (direction == Direction.Both || direction == Direction.Horizontal)
                cacheRows = PolyphaseCache.Build(M, this.m, this.window);

            if (direction == Direction.Both)
            {
                // 1) Left-multiply: U^H · A  → shape: (2N × M)
                var tmp = new Complex32[2 * N, M];

                Parallel.For(0, M, j =>
                {
                    var col = new Complex32[N];
                    for (int i = 0; i < N; i++) col[i] = A[i, j];

                    var tr = FWHT(col, cacheCols); // length 2N
                    for (int i = 0; i < 2 * N; i++) tmp[i, j] = tr[i];
                });

                // 2) Right-multiply by V via conjugation trick:
                //    row * V = conj( FWHT( conj(row) ) )  → shape: (2N × 2M)
                var B = new Complex32[2 * N, 2 * M];

                Parallel.For(0, 2 * N, i =>
                {
                    var row = new Complex32[M];
                    for (int j = 0; j < M; j++) row[j] = tmp[i, j];

                    // conj → FWHT → conj
                    for (int j = 0; j < M; j++) row[j] = new Complex32(row[j].Real, -row[j].Imag);
                    var tr = FWHT(row, cacheRows); // V^H * conj(row)

                    for (int j = 0; j < 2 * M; j++)
                    {
                        var c = tr[j];
                        B[i, j] = new Complex32(c.Real, -c.Imag); // conj(...)
                    }
                });

                return B;
            }

            if (direction == Direction.Vertical)
            {
                // B = U^H · A  → (2N × M)
                var B = new Complex32[2 * N, M];

                Parallel.For(0, M, j =>
                {
                    var col = new Complex32[N];
                    for (int i = 0; i < N; i++) col[i] = A[i, j];

                    var tr = FWHT(col, cacheCols); // 2N
                    for (int i = 0; i < 2 * N; i++) B[i, j] = tr[i];
                });

                return B;
            }
            else // Horizontal
            {
                // B = A · V  → (N × 2M)  (via the same conjugation trick per row)
                var B = new Complex32[N, 2 * M];

                Parallel.For(0, N, i =>
                {
                    var row = new Complex32[M];
                    for (int j = 0; j < M; j++) row[j] = A[i, j];

                    for (int j = 0; j < M; j++) row[j] = new Complex32(row[j].Real, -row[j].Imag);
                    var tr = FWHT(row, cacheRows); // V^H * conj(row)

                    for (int j = 0; j < 2 * M; j++)
                    {
                        var c = tr[j];
                        B[i, j] = new Complex32(c.Real, -c.Imag); // conj(...)
                    }
                });

                return B;
            }
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public override Complex32[,] Backward(Complex32[,] B)
        {
            int N2 = B.GetLength(0); // possibly 2N
            int M2 = B.GetLength(1); // possibly 2M

            PolyphaseCache cacheCols = null; // for target N = N2/2 (vertical inverse)
            PolyphaseCache cacheRows = null; // for target M = M2/2 (horizontal inverse)

            if (direction == Direction.Both || direction == Direction.Vertical)
                cacheCols = PolyphaseCache.Build(N2 / 2, this.m, this.window);

            if (direction == Direction.Both || direction == Direction.Horizontal)
                cacheRows = PolyphaseCache.Build(M2 / 2, this.m, this.window);

            if (direction == Direction.Both)
            {
                // Inverse order of Forward(Both):
                // 1) Undo right multiply (· V) per row using the conjugation trick:
                //    row = conj( IFWHT( conj(B_row) ) )  → shape: (N2 × M2/2)
                var tmp = new Complex32[N2, M2 / 2];

                Parallel.For(0, N2, i =>
                {
                    var row = new Complex32[M2];
                    for (int j = 0; j < M2; j++) row[j] = B[i, j];

                    for (int j = 0; j < M2; j++) row[j] = new Complex32(row[j].Real, -row[j].Imag); // conj(B_row)
                    var inv = IFWHT(row, cacheRows); // length M2/2

                    for (int j = 0; j < M2 / 2; j++)
                    {
                        var c = inv[j];
                        tmp[i, j] = new Complex32(c.Real, -c.Imag); // conj(...)
                    }
                });

                // 2) Undo left multiply (U^H ·) per column: 2N → N
                var A = new Complex32[N2 / 2, M2 / 2];

                Parallel.For(0, M2 / 2, j =>
                {
                    var col = new Complex32[N2];
                    for (int i = 0; i < N2; i++) col[i] = tmp[i, j];

                    var inv = IFWHT(col, cacheCols); // length N2/2
                    for (int i = 0; i < N2 / 2; i++) A[i, j] = inv[i];
                });

                return A;
            }

            if (direction == Direction.Vertical)
            {
                // A = (U^H)^{-1} · B : columns 2N → N
                var A = new Complex32[N2 / 2, M2];

                Parallel.For(0, M2, j =>
                {
                    var col = new Complex32[N2];
                    for (int i = 0; i < N2; i++) col[i] = B[i, j];

                    var inv = IFWHT(col, cacheCols); // N2/2
                    for (int i = 0; i < N2 / 2; i++) A[i, j] = inv[i];
                });

                return A;
            }
            else // Horizontal
            {
                // A = B · V^{-1} : rows 2M → M (via conjugation trick)
                var A = new Complex32[N2, M2 / 2];

                Parallel.For(0, N2, i =>
                {
                    var row = new Complex32[M2];
                    for (int j = 0; j < M2; j++) row[j] = B[i, j];

                    for (int j = 0; j < M2; j++) row[j] = new Complex32(row[j].Real, -row[j].Imag); // conj(B_row)
                    var inv = IFWHT(row, cacheRows); // M2/2

                    for (int j = 0; j < M2 / 2; j++)
                    {
                        var c = inv[j];
                        A[i, j] = new Complex32(c.Real, -c.Imag); // conj(...)
                    }
                });

                return A;
            }
        }

        #endregion

        #region Private voids

        /// <summary>
        /// Forward fast Weyl–Heisenberg transform (2-channel packing, length = 2N).
        /// 
        /// Input:
        ///   A ∈ ℂ^N,  N = M * L,  M is even.
        /// 
        /// Output:
        ///   B ∈ ℂ^{2N}, split as:
        ///     B_main[u] = B[u + 0],   u = l*M + k
        ///     B_half[u] = B[u + N],   u = l*M + k
        ///
        /// This matches the column-wise “slow” matrix reference G ∈ ℂ^{N×2N} with columns:
        ///   main: g[i] * exp( j*2π*k/M * (n - M/4) ) / √2
        ///   half: j*g[j] * exp( j*2π*k/M * (n - M/4) ) / √2
        ///
        /// Analysis (i.e., G^H A) realized via:
        ///   - polyphase split over residue n0 (mod M) and FFT_L along r,
        ///   - frequency-domain correlations with conj(S_hat/T_hat),
        ///   - carry-phase compensation for the half branch when (n0 + M/2) wraps,
        ///   - DFT_M over n0 (via forward FFT over n0) + quarter-phase e^{+jπk/2},
        ///   - normalization gain = √M / (√N · √2) = 1 / (√L · √2).
        ///
        /// Notes:
        ///   Caches C.S_hat, C.T_hat must be built from the SAME orthonormalized window
        ///   as used by the slow matrix builder, otherwise values will differ by window scaling.
        /// </summary>
        /// <param name="A">Input signal (length N)</param>
        /// <param name="C">Polyphase cache (built for N,M,window)</param>
        /// <returns>B of length 2N: main (0..N-1) and half (N..2N-1) branches</returns>
        internal static Complex32[] FWHT(Complex32[] A, PolyphaseCache C)
        {
            int N = A.Length;
            if (!Maths.IsPower(N, 2)) throw new Exception("Dimension of the signal must be a power of 2");

            var B = new Complex32[2 * N];

            int Mloc = C.M;
            if (!Maths.IsEven(Mloc)) throw new Exception("M must be even");

            int L = N / Mloc;
            if (L * Mloc != N) throw new Exception("N must be divisible by M");

            // Cached FFT_L spectra of window polyphase components (built from orthonormalized g):
            //   S_hat[n0,q] = FFT_L{ g[r*M + n0] }_r
            //   T_hat[n0,q] = FFT_L{ g[r*M + (n0+M/2) mod M] }_r
            // Dimensions: [M, L].
            var S_hat = C.S_hat;
            var T_hat = C.T_hat;

            // 1) Polyphase split over residue n0 (mod M). For each n0 we FFT along r (length L):
            //    Xhat[n0,q] = FFT_L{ A[r*M + n0] }_r
            var Xhat = new Complex32[Mloc, L];
            var tmp = new Complex32[L];

            for (int n0 = 0; n0 < Mloc; n0++)
            {
                for (int r = 0; r < L; r++)
                    tmp[r] = A[r * Mloc + n0];

                FFT(tmp, false); // forward FFT along r (no extra scaling)

                for (int q = 0; q < L; q++)
                    Xhat[n0, q] = tmp[q];
            }

            // 2) Frequency-domain correlations to accumulate over time-shift l:
            //
            // main branch:
            //   Cmain[n0,l] = IFFT_L{ conj(S_hat[n0,q]) * Xhat[n0,q] }_q
            //
            // half branch:
            //   Chalf[n0,l] = IFFT_L{ conj(T_hat[n0,q]) * Xhat[n0,q] * phase_carry(q) }_q
            //
            // carry phase:
            //   when (n0 + M/2) >= M, the half-branch index wraps and induces a +1 shift in r.
            //   A +1 shift in r corresponds in frequency to multiplication by exp(-j*2π*q/L).
            var Cmain = new Complex32[Mloc, L];
            var Chalf = new Complex32[Mloc, L];

            for (int n0 = 0; n0 < Mloc; n0++)
            {
                // main: correlate in frequency and go back to l-domain
                for (int q = 0; q < L; q++)
                    tmp[q] = S_hat[n0, q].Conjugate * Xhat[n0, q];

                FFT(tmp, true); // inverse FFT over r → correlation over l
                for (int l = 0; l < L; l++)
                    Cmain[n0, l] = tmp[l];

                // half: same, but with carry-phase if (n0 + M/2) wraps
                int carry = ((n0 + Mloc / 2) >= Mloc) ? 1 : 0;

                for (int q = 0; q < L; q++)
                {
                    float ang = 2f * Maths.Pi * q * carry / L;
                    var shiftPhase = Maths.Exp(-Complex32.I * ang); // exp(-j*2π*q/L) when carry=1

                    tmp[q] = T_hat[n0, q].Conjugate * Xhat[n0, q] * shiftPhase;
                }

                FFT(tmp, true); // inverse FFT over r
                for (int l = 0; l < L; l++)
                    Chalf[n0, l] = tmp[l];
            }

            // 3) Assemble over frequency shifts k for each time shift l.
            //
            // Over n0 we need a DFT with positive exponent e^{+j2π kn0/M}.
            // We realize it by a forward FFT over n0 (no extra *M here).
            // Then we apply the quarter-phase e^{+jπk/2} (accounts for the (n - M/4) shift).
            //
            // Final normalization:
            //   gain = √M / (√N · √2) = 1 / (√L · √2).
            //
            // Output packing (two channels):
            //   B_main[u] =  P * gain
            //   B_half[u] = (-j) Q * gain
            // where P, Q are the assembled complex contributions of main/half branches.
            var Sp_main = new Complex32[Mloc];
            var Sp_half = new Complex32[Mloc];

            float gain = Maths.Sqrt(Mloc) / (Maths.Sqrt(N) * Maths.Sqrt(2f)); // = 1/(√L·√2)

            for (int l = 0; l < L; l++)
            {
                // collect l-th rows across n0
                for (int n0 = 0; n0 < Mloc; n0++)
                {
                    Sp_main[n0] = Cmain[n0, l];
                    Sp_half[n0] = Chalf[n0, l];
                }

                // forward DFT over n0 via FFT (positive-exponent convention)
                FFT(Sp_main, false);
                FFT(Sp_half, false);

                for (int k = 0; k < Mloc; k++)
                {
                    var phase = PhasePlusPiOver2(k); // e^{+jπk/2}
                    var P = phase * Sp_main[k];      // main contribution
                    var Q = phase * Sp_half[k];      // half contribution

                    int u = l * Mloc + k;

                    // Two-channel output that matches the slow matrix reference (G^H A):
                    B[u + 0] = P * gain;                  // main channel
                    B[u + N] = -Complex32.I * Q * gain;   // half channel (−j factor)
                }
            }

            return B;
        }

        /// <summary>
        /// Inverse fast Weyl–Heisenberg transform (2-channel input, length = 2N).
        ///
        /// Input:
        ///   B ∈ ℂ^{2N}, split as:
        ///     B_main[u] = B[u + 0],   u = l*M + k
        ///     B_half[u] = B[u + N],   u = l*M + k
        ///
        /// This inverts the forward packing:
        ///   forward: B_main =  P * gain
        ///            B_half = (-j) Q * gain
        ///   where gain = √M / (√N · √2) = 1 / (√L · √2).
        ///
        /// Inverse mapping:
        ///   P = B_main / gain
        ///   Q =  j * B_half / gain
        ///
        /// Then we undo the quarter-phase e^{+jπk/2}, invert the DFT over n₀,
        /// apply the adjoint of the frequency-domain correlations, and finally
        /// invert the polyphase FFT over r to reconstruct A.
        /// </summary>
        /// <param name="B">Input coefficients of length 2N (main first, then half)</param>
        /// <param name="C">Polyphase cache (built for N, M, window)</param>
        /// <returns>Reconstructed signal A of length N</returns>
        internal static Complex32[] IFWHT(Complex32[] B, PolyphaseCache C)
        {
            int N = C.N;
            if (!Maths.IsPower(N, 2)) throw new Exception("Dimension of the signal must be a power of 2");
            if (B.Length != 2 * N) throw new Exception("Expect 2N coefficients (main + half)");

            int Mloc = C.M;
            if (!Maths.IsEven(Mloc)) throw new Exception("M must be even");

            int L = N / Mloc;
            if (L * Mloc != N) throw new Exception("N must be divisible by M");

            // Cached FFT_L spectra of window polyphase components:
            var S_hat = C.S_hat; // [M, L]
            var T_hat = C.T_hat; // [M, L]

            // ----- 1) Undo k-assembly (for each l) -----
            //
            // From B_main/B_half recover P,Q, remove quarter-phase, and
            // invert the DFT over n0 (IFFT_M) to obtain Cmain[:,l], Chalf[:,l].
            var Cmain = new Complex32[Mloc, L];
            var Chalf = new Complex32[Mloc, L];

            var Y_main = new Complex32[Mloc]; // spectra over k for a fixed l
            var Y_half = new Complex32[Mloc];

            // gain used in the forward:
            //   gain = 1 / (√L · √2)
            // hence inverse gain:
            float invGain = 1.0f / Maths.Sqrt(2 * L);

            for (int l = 0; l < L; l++)
            {
                // rebuild Y_main[k], Y_half[k]
                for (int k = 0; k < Mloc; k++)
                {
                    int u = l * Mloc + k;

                    var bMain = B[u + 0];
                    var bHalf = B[u + N];

                    // Invert forward packing:
                    //   P =  B_main / gain
                    //   Q =  j * B_half / gain
                    var P = bMain * invGain;
                    var Q = Complex32.I * bHalf * invGain;

                    // Remove quarter-phase e^{+jπk/2} → multiply by its conjugate
                    var phase = PhasePlusPiOver2(k);
                    var phaseConj = new Complex32(phase.Real, -phase.Imag);

                    Y_main[k] = phaseConj * P;
                    Y_half[k] = phaseConj * Q;
                }

                // Invert DFT over n0: IFFT_M.
                // Our FFT(true) applies an internal 1/√M; compensate by √M to get unit gain.
                FFT(Y_main, true);
                FFT(Y_half, true);

                float cM = Maths.Sqrt(Mloc); // compensates the internal 1/√M from IFFT
                for (int n0 = 0; n0 < Mloc; n0++)
                {
                    Cmain[n0, l] = Y_main[n0] * cM;
                    Chalf[n0, l] = Y_half[n0] * cM;
                }
            }

            // ----- 2) Adjoint of frequency-domain correlations (over l) -----
            //
            // Forward used:
            //   Cmain = IFFT_L{ conj(S_hat) * Xhat }
            //   Chalf = IFFT_L{ conj(T_hat) * Xhat * phase_carry }
            //
            // Adjoint (backward) therefore uses:
            //   Xhat += S_hat * FFT_L{ Cmain }
            //   Xhat += T_hat * FFT_L{ Chalf } * conj(phase_carry)
            var Xhat = new Complex32[Mloc, L];
            var bufL = new Complex32[L];
            var bufL2 = new Complex32[L];

            for (int n0 = 0; n0 < Mloc; n0++)
            {
                // FFT_L over l
                for (int l = 0; l < L; l++) bufL[l] = Cmain[n0, l];
                FFT(bufL, false);

                for (int l = 0; l < L; l++) bufL2[l] = Chalf[n0, l];
                FFT(bufL2, false);

                int carry = ((n0 + Mloc / 2) >= Mloc) ? 1 : 0;

                for (int q = 0; q < L; q++)
                {
                    // main branch adjoint
                    Xhat[n0, q] += S_hat[n0, q] * bufL[q];

                    // half branch adjoint (conjugate of forward’s carry-phase)
                    if (carry != 0)
                    {
                        float ang = 2f * Maths.Pi * q * carry / L;
                        var shiftPhaseAdj = Maths.Exp(Complex32.I * ang); // exp(+j*2π*q/L)
                        Xhat[n0, q] += T_hat[n0, q] * shiftPhaseAdj * bufL2[q];
                    }
                    else
                    {
                        Xhat[n0, q] += T_hat[n0, q] * bufL2[q];
                    }
                }
            }

            // ----- 3) Inverse polyphase recomposition (over r) -----
            //
            // For each residue n0, take IFFT_L over q to obtain samples along r:
            //   A[r*M + n0] = IFFT_L{ Xhat[n0, q] }_q
            var Arec = new Complex32[N];

            for (int n0 = 0; n0 < Mloc; n0++)
            {
                for (int q = 0; q < L; q++)
                    bufL[q] = Xhat[n0, q];

                FFT(bufL, true); // inverse FFT along r (library applies internal 1/√L)

                for (int r = 0; r < L; r++)
                    Arec[r * Mloc + n0] = bufL[r];
            }

            return Arec;
        }

        /// <summary>
        /// Defines the polyphase cache for the fast Weyl–Heisenberg transform.
        /// <para>
        /// This cache stores precomputed FFTs of the polyphase components of the
        /// orthogonalized analysis window. It is used in both forward and backward
        /// transforms to avoid recomputing window-related data for each call.
        /// </para>
        /// <para>
        /// Parameters:
        /// <list type="bullet">
        /// <item><description><see cref="N"/> – total signal length</description></item>
        /// <item><description><see cref="M"/> – number of frequency shifts (must be even)</description></item>
        /// <item><description><see cref="L"/> – number of time shifts, L = N / M</description></item>
        /// <item><description><see cref="S_hat"/> – FFT<sub>L</sub> of the main polyphase branch
        /// s<sub>n0</sub>[r] = g[r*M + n0], dimensions [M, L]</description></item>
        /// <item><description><see cref="T_hat"/> – FFT<sub>L</sub> of the half-shifted polyphase branch
        /// t<sub>n0</sub>[r] = g[r*M + (n0 + M/2) mod M], dimensions [M, L]</description></item>
        /// </list>
        /// </para>
        /// <remarks>
        /// The cache is specific to a given combination of (N, M, window function).
        /// If the window changes, the cache must be rebuilt.
        /// </remarks>
        /// </summary>
        internal sealed class PolyphaseCache
        {
            /// <summary>
            /// Total signal length.
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

            /// <summary>
            /// FFT<sub>L</sub> of the main branch polyphase components [M, L].
            /// </summary>
            public readonly Complex32[,] S_hat;

            /// <summary>
            /// FFT<sub>L</sub> of the half-shifted branch polyphase components [M, L].
            /// </summary>
            public readonly Complex32[,] T_hat;

            /// <summary>
            /// Creates a new instance of the polyphase cache.
            /// </summary>
            /// <param name="N">Total signal length</param>
            /// <param name="M">Number of frequency shifts</param>
            /// <param name="L">Number of time shifts</param>
            /// <param name="S">FFT of the main branch polyphase components</param>
            /// <param name="T">FFT of the half-shifted branch polyphase components</param>
            public PolyphaseCache(int N, int M, int L, Complex32[,] S, Complex32[,] T)
            {
                this.N = N;
                this.M = M;
                this.L = L;
                this.S_hat = S;
                this.T_hat = T;
            }

            /// <summary>
            /// Precomputes FFT_L of the polyphase components of the orthogonalized window.
            /// These cached arrays S_hat and T_hat are used in both Forward() and Backward()
            /// to avoid recomputing window FFTs for every transform call.
            /// </summary>
            /// <param name="N">Total signal length</param>
            /// <param name="Mloc">Number of frequency shifts M (must be even)</param>
            /// <param name="window">Windows function</param>
            /// <returns>Polyphase cache</returns>
            public static PolyphaseCache Build(int N, int Mloc, IWindow window)
            {
                // If cache for the given N is already computed, reuse it
                if (!Maths.IsPower(N, 2)) throw new Exception("Dimension of the signal must be a power of 2");
                if (!Maths.IsEven(Mloc)) throw new Exception("M must be even");
                int L = N / Mloc;
                if (L * Mloc != N) throw new Exception("N must be divisible by M");

                // Step 1: Generate the initial WH analysis window g0 of length N
                //         from the current IWindow object (e.g., Gaussian, Hanning, etc.)
                var g0 = FastWeylHeisenbergTransform.GetPacket(window, N);

                // Step 2: Orthogonalize g0 using Zak-domain orthogonalization
                //         to produce a WH-orthonormal window g.
                //         This matches the behavior of Matrix(..., true) in the slow implementation.
                var zakOrth = new FastZakTransform(Mloc);
                var g = zakOrth.Orthogonalize(g0); // real-valued array of length N

                // Allocate caches for the FFTs of polyphase components:
                // S_hat[n0, q] = FFT_L over r of polyphase component g[r*M + n0]
                // T_hat[n0, q] = FFT_L over r of polyphase component g[r*M + (n0 + M/2) % M]
                var S_hat = new Complex32[Mloc, L];
                var T_hat = new Complex32[Mloc, L];

                var s = new Complex32[L]; // temporary buffer for length-L FFT

                for (int n0 = 0; n0 < Mloc; n0++)
                {
                    // --- Main branch polyphase component ---
                    // s[r] = g[r*M + n0], r = 0..L-1
                    for (int r = 0; r < L; r++)
                        s[r] = new Complex32(g[r * Mloc + n0], 0);

                    // Forward FFT along r (length L)
                    FFT(s, false);

                    // Store as S_hat[n0, :]
                    for (int q = 0; q < L; q++)
                        S_hat[n0, q] = s[q];

                    // --- Half-shifted branch polyphase component ---
                    // Index n1 is (n0 + M/2) mod M — frequency index shifted by half the band.
                    int n1 = (n0 + Mloc / 2) % Mloc;

                    for (int r = 0; r < L; r++)
                        s[r] = new Complex32(g[r * Mloc + n1], 0);

                    // Forward FFT along r
                    FFT(s, false);

                    // Store as T_hat[n0, :]
                    for (int q = 0; q < L; q++)
                        T_hat[n0, q] = s[q];
                }

                // return cache
                return new PolyphaseCache(N, Mloc, L, S_hat, T_hat);
            }
        }

        /// <summary>
        /// Returns the complex phase factor e^{+j * π * k / 2}.
        /// This term is 4-periodic in k and takes only four distinct values: 1, +j, −1, −j.
        /// 
        /// Context:
        /// • Quarter-period phase factors appear in WH transforms, DFT twiddle factors,
        ///   and modulation/demodulation stages.
        /// • In the Weyl–Heisenberg transform, this multiplier is used during n₀→k assembly
        ///   with a positive-exponent DFT to align phases with the matrix reference.
        /// 
        /// Relationship:
        ///   PhasePlusPiOver2(k) = conj(PhaseMinusPiOver2(k)).
        /// 
        /// Mapping (k mod 4):
        ///   0 →  1   (  0°)
        ///   1 → +j   (+90°)
        ///   2 → −1   (180°)
        ///   3 → −j   (−90°)
        /// </summary>
        /// <param name="k">Frequency index (integer).</param>
        /// <returns>Complex value of e^{+j * π * k / 2}.</returns>
        private static Complex32 PhasePlusPiOver2(int k)
        {
            // Fast k % 4 using bitwise AND with 3 (0b11).
            return (k & 3) switch
            {
                0 => new Complex32(+1, 0),  //   0°
                1 => new Complex32(0, +1),  // +90°
                2 => new Complex32(-1, 0),  // 180°
                _ => new Complex32(0, -1),  // −90°
            };
        }

        /// <summary>
        /// Fast Fourier transform.
        /// </summary>
        /// <param name="a">Input</param>
        /// <param name="inverse">Iverse or not</param>
        private static void FFT(Complex32[] a, bool inverse)
        {
            var n = a.Length;

            if (inverse)
            {
                // Perform the inverse FFT (frequency → time domain)
                var b = fastFourierTransform.Backward(a);

                // Apply orthonormal scaling by sqrt(1 / N) to match analysis/synthesis norms
                float inv = Maths.Sqrt(1f / n);

                for (int i = 0; i < n; i++)
                {
                    a[i] = b[i] * inv;
                }
            }
            else
            {
                // Perform the forward FFT (time → frequency domain)
                var b = fastFourierTransform.Forward(a);

                for (int i = 0; i < n; i++)
                {
                    a[i] = b[i];
                }
            }
        }

        /// <summary>
        /// UMapx fast Fourier transform.
        /// </summary>
        private static readonly FastFourierTransform fastFourierTransform = new FastFourierTransform(false, Direction.Vertical);

        #endregion
    }
}
