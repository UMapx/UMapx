using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines fast Weyl-Heisenberg transform.
    /// </summary>
    /// <remarks>
    /// The class represents a computationally efficient implementation of one-dimensional and two-dimensional discrete orthogonal
    /// Weyl-Heisenberg transforms. This implementation was designed and developed by Valery Asiryan, Yerevan, Armenia (2025).
    /// </remarks>
    public partial class FastWeylHeisenbergTransform
    {
        #region Internal methods and classes

        /// <summary>
        /// Forward fast Weyl–Heisenberg transform (2-channel packing, output length = 2N).
        /// 
        /// Signal model and factorization:
        ///   • N = M · L,  M is even. Index the signal as n = r·M + n0 with residues n0∈[0..M-1], r∈[0..L-1].
        ///   • Polyphase split over n0 and length-L FFT over r:
        ///       Xhat[n0,q] = FFT_L{ A[r*M + n0] }_r.
        ///   • Frequency-domain correlations (analysis) with orthonormalized window spectra:
        ///       main : C_main[n0,l] = IFFT_L{ conj(S_hat[n0,q]) * Xhat[n0,q] }_q,
        ///       half : C_half[n0,l] = IFFT_L{ conj(T_hat[n0,q]) * Xhat[n0,q] * φ_carry(q) }_q,
        ///     where φ_carry(q) = exp(−j·2π q/L) iff (n0 + M/2) wraps modulo M (time shift +1 in r).
        ///   • Across residues n0, assemble frequency shifts k by a forward DFT_M with positive exponent
        ///     and apply the quarter-period phase e^{+jπk/2} to account for the (n − M/4) centering.
        ///   • Two-channel packing:
        ///       B_main[u] =  P(l,k) · gain,
        ///       B_half[u] = −j · Q(l,k) · gain,
        ///     with u = l·M + k and gain = √M / √N = 1/√L (orthonormal convention).
        ///
        /// Exactness:
        ///   Matches the column-wise “slow” reference matrix G ∈ ℂ^{N×2N} built from the SAME
        ///   Zak-orthonormalized analysis window g used to compute S_hat/T_hat in PolyphaseCache.
        /// 
        /// Complexity:
        ///   O(M·L·log L) for the r-axis FFT/IFFT correlations + O(L·M·log M) for the n0→k assembly.
        ///   Memory uses O(L·M) only for two row-major workspaces C_main and C_half (no Xhat buffer).
        /// </summary>
        /// <param name="A">Input signal A ∈ ℂ^N, N = M·L</param>
        /// <param name="C">Polyphase cache holding S_hat and T_hat computed from the orthonormal window</param>
        /// <returns>B ∈ ℂ^{2N}: main in B[0..N-1], half in B[N..2N-1], packed as u = l·M + k</returns>
        internal static Complex32[] FWHT(Complex32[] A, PolyphaseCache C)
        {
            int N = A.Length;
            int M = C.M;
            int L = N / M;
            if (L * M != N) throw new ArgumentException("N must be divisible by M");

            var S_hat = C.S_hat; // [M][L]  spectra of g[r*M + n0] over r
            var T_hat = C.T_hat; // [M][L]  spectra of g[r*M + (n0 + M/2) mod M] over r

            // Precompute carry-phase φ_carry(q) = exp(−j·2π q/L).
            // It models a +1 cyclic shift in r (time) caused by index wrap when (n0 + M/2) ≥ M.
            var shiftL = new Complex32[L];
            for (int q = 0; q < L; q++)
            {
                float ang = 2f * Maths.Pi * q / L;
                shiftL[q] = Maths.Exp(-Complex32.I * ang);
            }

            // Row-major accumulation buffers:
            //   Cmain_rows[l][n0]  holds C_main[n0,l],
            //   Chalf_rows[l][n0]  holds C_half[n0,l].
            // We store by rows to run FFT_M directly in-place per row (no extra copies).
            var Cmain_rows = new Complex32[L][];
            var Chalf_rows = new Complex32[L][];
            for (int l = 0; l < L; l++)
            {
                Cmain_rows[l] = new Complex32[M];
                Chalf_rows[l] = new Complex32[M];
            }

            // Length-L temporaries reused for all residues.
            var X = new Complex32[L]; // X[q]  = FFT_L polyphase spectrum for residue n0
            var Y = new Complex32[L]; // work  = pointwise products in frequency

            // ----- Stage 1: polyphase → FFT_L over r → frequency correlations → IFFT_L to l domain -----
            for (int n0 = 0; n0 < M; n0++)
            {
                // Polyphase gather along r for this residue n0.
                for (int r = 0; r < L; r++)
                    X[r] = A[r * M + n0];

                // X[q] ← FFT_L{ A[r*M + n0] }_r  (forward FFT, no extra scaling).
                FFT(X, false);

                // MAIN branch: Y[q] = conj(S_hat[n0,q]) * X[q].
                var Sh = S_hat[n0];
                for (int q = 0; q < L; q++)
                    Y[q] = Sh[q].Conjugate * X[q];

                // IFFT_L over q: gives correlation sequence over time shifts l.
                FFT(Y, true); // inverse is orthonormal (applies 1/√L internally)
                for (int l = 0; l < L; l++)
                    Cmain_rows[l][n0] = Y[l];

                // HALF branch:
                // If (n0 + M/2) wraps past M, that implies a +1 shift in r, i.e. multiply by exp(−j·2πq/L).
                var Th = T_hat[n0];
                bool carry = ((n0 + (M >> 1)) >= M);
                if (carry)
                {
                    for (int q = 0; q < L; q++)
                        Y[q] = Th[q].Conjugate * X[q] * shiftL[q];
                }
                else
                {
                    for (int q = 0; q < L; q++)
                        Y[q] = Th[q].Conjugate * X[q];
                }

                // IFFT_L → half branch correlations over l.
                FFT(Y, true);
                for (int l = 0; l < L; l++)
                    Chalf_rows[l][n0] = Y[l];
            }

            // ----- Stage 2: per-row assembly over k (DFT_M with + exponent) + quarter-phase + packing -----
            var B = new Complex32[2 * N];
            float gain = Maths.Sqrt(M) / Maths.Sqrt(N); // = 1/√L (unitary end-to-end)

            for (int l = 0; l < L; l++)
            {
                var rowM = Cmain_rows[l]; // length M
                var rowH = Chalf_rows[l]; // length M

                // Forward DFT over residues n0 (positive exponent) via FFT.
                FFT(rowM, false);
                FFT(rowH, false);

                // Apply quarter-phase e^{+jπk/2} (4-periodic: 1, +j, −1, −j) and pack outputs.
                for (int k = 0; k < M; k++)
                {
                    var phase = PhasePlusPiOver2(k); // e^{+jπk/2}

                    // P = e^{+jπk/2} · FFT_M{C_main}[k]
                    // Q = e^{+jπk/2} · FFT_M{C_half}[k]
                    var P = phase * rowM[k];
                    var Q = phase * rowH[k];

                    int u = l * M + k;

                    // Two-channel packing consistent with the slow matrix reference:
                    //   main:  B[u]     =  P · gain,
                    //   half:  B[u + N] = (−j) · Q · gain.
                    B[u] = P * gain;
                    B[u + N] = -Complex32.I * Q * gain;
                }
            }

            return B;
        }

        /// <summary>
        /// Inverse fast Weyl–Heisenberg transform (synthesis), inverting the forward packing exactly.
        /// 
        /// Given B_main and B_half (packed by u = l·M + k):
        ///   • Undo −j on the half branch and the quarter-phase e^{+jπk/2}, recover the spectra
        ///     over k, then IFFT_M to get the per-residue correlation sequences C_main[:,l], C_half[:,l].
        ///   • FFT_L over l to return to the r-frequency axis (q), then apply the adjoint frequency
        ///     correlations (synthesis):  Xhat[n0,q] = S_hat[n0,q]·FFT{C_main} + T_hat[n0,q]·FFT{C_half}·conj(φ_carry).
        ///   • IFFT_L over q and interleave residues to reconstruct A[n] = A[r·M + n0].
        ///
        /// Normalization tracks the forward:
        ///   forward gain = √M/√N = 1/√L and inverse uses the corresponding 1/gain as part of the unpacking
        ///   (here kept as <c>invGain = 1/(2√L)</c> to match the user’s established convention and matrix reference).
        ///
        /// Complexity:
        ///   O(L·M·log M) for IFFT_M over k + O(M·L·log L) for synthesis along the r-axis.
        /// </summary>
        /// <param name="B">Input coefficients B ∈ ℂ^{2N}: main in [0..N-1], half in [N..2N-1]</param>
        /// <param name="C">Polyphase cache with the same orthonormal window spectra S_hat/T_hat</param>
        /// <returns>Reconstructed signal A ∈ ℂ^N</returns>
        internal static Complex32[] IFWHT(Complex32[] B, PolyphaseCache C)
        {
            int N = C.N;
            int M = C.M;
            int L = N / M;

            if (B.Length != 2 * N) throw new ArgumentException("Expect 2N coefficients");
            if (L * M != N) throw new ArgumentException("N must be divisible by M");

            var S_hat = C.S_hat; // [M][L]
            var T_hat = C.T_hat; // [M][L]

            // Conjugate of the forward carry-phase:
            //   conj( exp(−j·2π q/L) ) = exp(+j·2π q/L).
            // Used when the half-branch residue wrapped in the forward.
            var shiftL_adj = new Complex32[L];
            for (int q = 0; q < L; q++)
            {
                float ang = 2f * Maths.Pi * q / L;
                shiftL_adj[q] = Maths.Exp(Complex32.I * ang);
            }

            // Row-major workspaces to hold C_main[:,l] and C_half[:,l] after undoing k-assembly.
            var Cmain_rows = new Complex32[L][];
            var Chalf_rows = new Complex32[L][];
            for (int l = 0; l < L; l++)
            {
                Cmain_rows[l] = new Complex32[M];
                Chalf_rows[l] = new Complex32[M];
            }

            var Y_main = new Complex32[M]; // spectra over k for a fixed l (main)
            var Y_half = new Complex32[M]; // spectra over k for a fixed l (half)

            // In the original codebase, the inverse unpacking used invGain = 1/(2√L).
            // We keep this constant to fully match the established normalization and the slow matrix.
            float invGain = 1.0f / (2 * Maths.Sqrt(L));

            // Compensate the internal 1/√M applied by IFFT_M in FFT(..., true) to keep unit gain.
            float cM = Maths.Sqrt(M);

            // ----- Stage 1: undo packing per (l,k) → IFFT_M over k to get row signals over n0 -----
            for (int l = 0; l < L; l++)
            {
                // For each frequency bin k, remove quarter-phase and −j on half branch,
                // and divide by the forward gain (embedded here via invGain).
                for (int k = 0; k < M; k++)
                {
                    int u = l * M + k;

                    var bMain = B[u];     // P · gain
                    var bHalf = B[u + N]; // (−j) · Q · gain

                    // Recover P and Q (inverse of forward packing):
                    //   P =      B_main / gain,
                    //   Q =  j · B_half / gain   (since forward multiplies half by −j).
                    var P = bMain * invGain;
                    var Q = Complex32.I * bHalf * invGain;

                    // Remove e^{+jπk/2} by multiplying with its conjugate.
                    var phase = PhasePlusPiOver2(k);
                    var phaseConj = new Complex32(phase.Real, -phase.Imag);

                    Y_main[k] = phaseConj * P;
                    Y_half[k] = phaseConj * Q;
                }

                // IFFT_M over k → row signals indexed by residue n0.
                FFT(Y_main, true);
                FFT(Y_half, true);

                // Compensate the internal 1/√M scaling to produce unit-gain rows.
                for (int n0 = 0; n0 < M; n0++)
                {
                    Cmain_rows[l][n0] = Y_main[n0] * cM;
                    Chalf_rows[l][n0] = Y_half[n0] * cM;
                }
            }

            // ----- Stage 2: synthesis along the r-axis (adjoint of analysis correlations) -----
            // For each residue n0:
            //   • take FFT_L over l for both main/half rows,
            //   • combine with S_hat/T_hat and conj(carry) in frequency (q),
            //   • IFFT_L over q to go back to r and interleave.
            var Arec = new Complex32[N];
            var bufMain = new Complex32[L];
            var bufHalf = new Complex32[L];
            var bufX = new Complex32[L];

            for (int n0 = 0; n0 < M; n0++)
            {
                // Gather column n0 across all l into contiguous buffers.
                for (int l = 0; l < L; l++)
                {
                    bufMain[l] = Cmain_rows[l][n0];
                    bufHalf[l] = Chalf_rows[l][n0];
                }

                // FFT_L over l → frequency domain along q.
                FFT(bufMain, false);
                FFT(bufHalf, false);

                // Frequency-domain synthesis:
                //   Xhat[n0,q] = S_hat[n0,q]*bufMain[q] + T_hat[n0,q]*bufHalf[q]*conj(φ_carry),
                //   where conj(φ_carry) = exp(+j·2π q/L) iff the forward half-branch wrapped.
                var Sh = S_hat[n0];
                var Th = T_hat[n0];
                bool carry = ((n0 + (M >> 1)) >= M);

                if (carry)
                {
                    for (int q = 0; q < L; q++)
                        bufX[q] = Sh[q] * bufMain[q] + Th[q] * shiftL_adj[q] * bufHalf[q];
                }
                else
                {
                    for (int q = 0; q < L; q++)
                        bufX[q] = Sh[q] * bufMain[q] + Th[q] * bufHalf[q];
                }

                // IFFT_L over q → r-domain samples for this residue n0, then de-interleave.
                FFT(bufX, true);
                for (int r = 0; r < L; r++)
                    Arec[r * M + n0] = bufX[r];
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
        /// </summary>
        /// <remarks>
        /// The cache is specific to a given combination of (N, M, window function).
        /// If the window changes, the cache must be rebuilt.
        /// </remarks>
        internal sealed class PolyphaseCache
        {
            /// <summary>
            /// Total signal length. Must satisfy N = M x L.
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
            /// FFT<sub>L</sub> of the main branch polyphase components [M][L].
            /// </summary>
            public readonly Complex32[][] S_hat;

            /// <summary>
            /// FFT<sub>L</sub> of the half-shifted branch polyphase components [M][L].
            /// </summary>
            public readonly Complex32[][] T_hat;

            /// <summary>
            /// Creates a new instance of the polyphase cache.
            /// </summary>
            /// <param name="N">Total signal length</param>
            /// <param name="M">Number of frequency shifts</param>
            /// <param name="L">Number of time shifts</param>
            /// <param name="S">FFT of the main branch polyphase components</param>
            /// <param name="T">FFT of the half-shifted branch polyphase components</param>
            public PolyphaseCache(int N, int M, int L, Complex32[][] S, Complex32[][] T)
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
                int L = N / Mloc;
                if (L * Mloc != N) throw new ArgumentException("N must be divisible by M");

                // Step 1: Generate the initial WH analysis window g0 of length N
                //         from the current IWindow object (e.g., Gaussian, Hanning, etc.)
                var g0 = FastWeylHeisenbergTransform.Packet(window, N);

                // Step 2: Orthogonalize g0 using Zak-domain orthogonalization
                //         to produce a WH-orthonormal window g.
                //         This matches the behavior of Matrix(..., true) in the slow implementation.
                var zakOrth = new FastZakTransform(Mloc);
                var g = zakOrth.Orthogonalize(g0); // real-valued array of length N

                // Allocate caches for the FFTs of polyphase components:
                // S_hat[n0, q] = FFT_L over r of polyphase component g[r*M + n0]
                // T_hat[n0, q] = FFT_L over r of polyphase component g[r*M + (n0 + M/2) % M]
                var S_hat = new Complex32[Mloc][];
                var T_hat = new Complex32[Mloc][];

                var s = new Complex32[L]; // temporary buffer for length-L FFT

                for (int n0 = 0; n0 < Mloc; n0++)
                {
                    // Initialize caches
                    S_hat[n0] = new Complex32[L];
                    T_hat[n0] = new Complex32[L];

                    // --- Main branch polyphase component ---
                    // s[r] = g[r*M + n0], r = 0..L-1
                    for (int r = 0; r < L; r++)
                        s[r] = new Complex32(g[r * Mloc + n0], 0);

                    // Forward FFT along r (length L)
                    FFT(s, false);

                    // Store as S_hat[n0, :]
                    for (int q = 0; q < L; q++)
                        S_hat[n0][q] = s[q];

                    // --- Half-shifted branch polyphase component ---
                    // Index n1 is (n0 + M/2) mod M — frequency index shifted by half the band.
                    int n1 = (n0 + Mloc / 2) % Mloc;

                    for (int r = 0; r < L; r++)
                        s[r] = new Complex32(g[r * Mloc + n1], 0);

                    // Forward FFT along r
                    FFT(s, false);

                    // Store as T_hat[n0, :]
                    for (int q = 0; q < L; q++)
                        T_hat[n0][q] = s[q];
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
        /// <param name="k">Frequency index (integer)</param>
        /// <returns>Complex value of e^{+j * π * k / 2}</returns>
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
        /// <param name="inverse">Inverse or not</param>
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
