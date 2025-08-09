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
        #region Private data
        private bool complex;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes fast Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N]</param>
        /// <param name="direction">Processing direction</param>
        public FastWeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical) : base(window, m, direction) { }
        /// <summary>
        /// Use the complex transform or not.
        /// <remarks>
        /// The algorithm will be fast only for the forward transform.
        /// </remarks>
        /// </summary>
        public bool Complex
        {
            get
            { 
                return this.complex; 
            }
            set
            { 
                this.complex = value; 
            }
        }
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
            return FWHT(A, cache, this.complex);
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            if (this.complex)
                return base.Backward(B);

            int N = B.Length;
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
            Complex32[,] B = (Complex32[,])A.Clone();
            int N = B.GetLength(0), M = B.GetLength(1);

            PolyphaseCache cacheCols = null;
            PolyphaseCache cacheRows = null;

            if (direction == Direction.Both || direction == Direction.Vertical)
                cacheCols = PolyphaseCache.Build(N, this.m, this.window);

            if (direction == Direction.Both || direction == Direction.Horizontal)
                cacheRows = PolyphaseCache.Build(M, this.m, this.window);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = FWHT(row, cacheRows, this.complex);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = FWHT(col, cacheCols, this.complex);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = FWHT(col, cacheCols, this.complex);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = FWHT(row, cacheRows, this.complex);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public override Complex32[,] Backward(Complex32[,] B)
        {
            if (this.complex)
                return base.Backward(B);

            Complex32[,] A = (Complex32[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            PolyphaseCache cacheCols = null;
            PolyphaseCache cacheRows = null;

            if (direction == Direction.Both || direction == Direction.Vertical)
                cacheCols = PolyphaseCache.Build(N, this.m, this.window);

            if (direction == Direction.Both || direction == Direction.Horizontal)
                cacheRows = PolyphaseCache.Build(M, this.m, this.window);

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = IFWHT(col, cacheCols);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = IFWHT(row, cacheRows);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = IFWHT(col, cacheCols);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = IFWHT(row, cacheRows);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        #endregion

        #region Private voids

        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <param name="C">Polyphase cache</param>
        /// <param name="complex">Complex or not</param>
        /// <returns>Array</returns>
        internal static Complex32[] FWHT(Complex32[] A, PolyphaseCache C, bool complex = false)
        {
            int N = A.Length;
            var B = new Complex32[N];

            // Use the complex transform or not
            // (works for complex signals)
            if (complex)
            {
                var Ar = new Complex32[N];
                var Ai = new Complex32[N];

                // Separate
                for (int n = 0; n < N; n++)
                {
                    Ar[n] = new Complex32(A[n].Real, 0);
                    Ai[n] = new Complex32(A[n].Imag, 0);
                }

                // Run the “real-input” fast WH twice
                var Br = FWHT(Ar, C);  // already matches matrix for real inputs
                var Bi = FWHT(Ai, C);

                // Combine: U^H(Ar + i Ai) = U^H Ar + i U^H Ai
                for (int u = 0; u < N; u++)
                {
                    B[u] = Br[u] + Complex32.I * Bi[u];
                }

                return B;
            }

            int Mloc = C.M;
            if (!Maths.IsEven(Mloc)) throw new Exception("M must be even");
            int L = N / Mloc;
            if (L * Mloc != N) throw new Exception("N must be divisible by M");

            // Precompute Zak/FFT caches of the (orthogonalized) window:
            // S_hat[n0, q] = FFT_L { g[r*M + n0] } over r
            // T_hat[n0, q] = FFT_L { g[r*M + (n0+M/2)%M] } over r
            // Dimensions: [M, L]
            var S_hat = C.S_hat; var T_hat = C.T_hat;

            // 1) Polyphase split of A by residue n0 (mod M), then FFT along r (length L):
            //    Xhat[n0, q] = FFT_L { A[r*M + n0] }_r
            var Xhat = new Complex32[Mloc, L];
            var tmp = new Complex32[L];

            for (int n0 = 0; n0 < Mloc; n0++)
            {
                for (int r = 0; r < L; r++)
                    tmp[r] = A[r * Mloc + n0];

                FFT(tmp, false); // forward FFT along r (no scaling assumed)

                for (int q = 0; q < L; q++)
                    Xhat[n0, q] = tmp[q];
            }

            // 2) Correlations along the time-shift index l via frequency domain:
            //    Cmain[n0, l]  = IFFT_L { conj(S_hat[n0,q]) * Xhat[n0,q] }_q
            //    Chalf[n0, l]  = IFFT_L { conj(T_hat[n0,q]) * Xhat[n0,q] * phase_carry(q) }_q
            //
            //    Here phase_carry(q) accounts for the wrap-around of the half-shifted branch:
            //    when (n0 + M/2) overflows modulo M, it induces a +1 shift in r, which is
            //    a multiplicative phase factor in frequency: exp(-j*2π*q/L) for carry=1.
            var Cmain = new Complex32[Mloc, L];
            var Chalf = new Complex32[Mloc, L];

            for (int n0 = 0; n0 < Mloc; n0++)
            {
                // --- main branch correlation ---
                for (int q = 0; q < L; q++)
                    tmp[q] = S_hat[n0, q].Conjugate * Xhat[n0, q];

                FFT(tmp, true); // inverse FFT along r (produces correlation sequence over l)

                for (int l = 0; l < L; l++)
                    Cmain[n0, l] = tmp[l];

                // --- half branch correlation with r-carry compensation ---
                // carry = 1 if (n0 + M/2) >= M, i.e., the half-shift in n0 crosses a block boundary in r
                int carry = ((n0 + Mloc / 2) >= Mloc) ? 1 : 0;

                for (int q = 0; q < L; q++)
                {
                    // Phase factor from a +carry circular shift in r-domain:
                    // shiftPhase(q) = exp(-j * 2π * q * carry / L)
                    float ang = 2f * Maths.Pi * q * carry / L;
                    var shiftPhase = Maths.Exp(-Complex32.I * ang);

                    // Correlation in frequency domain for half branch (note conj on T_hat):
                    tmp[q] = T_hat[n0, q].Conjugate * Xhat[n0, q] * shiftPhase;
                }

                FFT(tmp, true); // inverse FFT along r

                for (int l = 0; l < L; l++)
                    Chalf[n0, l] = tmp[l];
            }

            // 3) Assemble over frequency shifts k for each time shift l:
            //    Over n0 we need a DFT with positive exponent => use IFFT_M (which yields +j sign)
            //    and then compensate its internal 1/M scaling by multiplying by M afterwards.
            //
            //    Finally apply the global phase exp(-j*pi*k/2) that matches the original basis
            //    (accounts for the (n - M/4) phase in the matrix version).
            var rowMain = new Complex32[Mloc];
            var rowHalf = new Complex32[Mloc];
            var Sp_main = new Complex32[Mloc];
            var Sp_half = new Complex32[Mloc];

            for (int l = 0; l < L; l++)
            {
                // Collect the l-th “row” across n0
                for (int n0 = 0; n0 < Mloc; n0++)
                {
                    rowMain[n0] = Cmain[n0, l];
                    rowHalf[n0] = Chalf[n0, l];
                }

                Array.Copy(rowMain, Sp_main, Mloc);
                Array.Copy(rowHalf, Sp_half, Mloc);

                // IFFT_M over n0 (positive-exponent DFT), then undo the internal /M by *M
                FFT(Sp_main, true);
                FFT(Sp_half, true);

                for (int k = 0; k < Mloc; k++)
                {
                    Sp_main[k] *= Mloc;
                    Sp_half[k] *= Mloc;
                }

                // Per-k phase alignment and channel packing (cos -> Re, sin -> Im)
                for (int k = 0; k < Mloc; k++)
                {
                    var phase = PhaseMinusPiOver2(k); // exp(-j*pi*k/2)

                    var P = phase * Sp_main[k]; // “cosine” branch contribution
                    var Q = phase * Sp_half[k]; // “sine”   branch contribution

                    // Match the matrix API scaling: divide by sqrt(N)
                    B[l * Mloc + k] = new Complex32(
                        P.Real / Maths.Sqrt(N),
                        Q.Imag / Maths.Sqrt(N)
                    );
                }
            }

            return B;
        }

        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <param name="C">Polyphase cache</param>
        /// <returns>Array</returns>
        internal static Complex32[] IFWHT(Complex32[] B, PolyphaseCache C)
        {
            int N = B.Length;
            int Mloc = C.M;
            if (!Maths.IsEven(Mloc)) throw new Exception("M must be even");
            int L = N / Mloc;
            if (L * Mloc != N) throw new Exception("N must be divisible by M");

            // Use the same window caches as in Forward():
            // S_hat/T_hat are FFT_L of the window’s polyphase components.
            var S_hat = C.S_hat; var T_hat = C.T_hat;

            // 1) Undo the k-assembly for each l:
            //    From B[l,k], rebuild Sp_main[k], Sp_half[k]:
            //      Forward had:  P = phase * Sp_main[k];  B.Re = Re(P)/sqrt(N)
            //                    Q = phase * Sp_half[k];  B.Im = Im(Q)/sqrt(N)
            //    Here we invert that mapping and remove the per-k phase.
            var Cmain = new Complex32[Mloc, L];
            var Chalf = new Complex32[Mloc, L];

            var Sp_main = new Complex32[Mloc];
            var Sp_half = new Complex32[Mloc];

            float s = Maths.Sqrt(N);

            for (int l = 0; l < L; l++)
            {
                for (int k = 0; k < Mloc; k++)
                {
                    var phase = PhaseMinusPiOver2(k);
                    var phaseConj = new Complex32(phase.Real, -phase.Imag); // conj(phase)

                    var b = B[l * Mloc + k];

                    // Invert packing: P := (b.Re * sqrt(N)) + 0i;   Q := i * (b.Im * sqrt(N))
                    var P = new Complex32(b.Real * s, 0);
                    var Q = new Complex32(0, b.Imag * s);

                    // Remove the per-k phase to obtain Sp_*[k]
                    Sp_main[k] = phaseConj * P;
                    Sp_half[k] = phaseConj * Q;
                }

                // Forward used: Sp_* = IFFT_M(row) * M
                // Therefore, row = FFT_M(Sp_*) / M  (no extra scaling terms)
                FFT(Sp_main, false); // FFT over n0 (length M)
                FFT(Sp_half, false);

                for (int n0 = 0; n0 < Mloc; n0++)
                {
                    Cmain[n0, l] = Sp_main[n0] * (1f / Maths.Sqrt(L * N));
                    Chalf[n0, l] = Sp_half[n0] * (1f / Maths.Sqrt(L * N));
                }
            }

            // 2) Undo the correlation step along l (apply adjoint operators in frequency):
            //    Forward: Cmain = IFFT_L{ conj(S_hat)*Xhat }, Chalf = IFFT_L{ conj(T_hat)*Xhat * phase_carry }
            //    Backward (adjoint): Xhat += S_hat * FFT_L{Cmain}, and
            //                        Xhat += T_hat * FFT_L{Chalf} * phase_carry_adj
            //
            //    For the half branch, the carry phase in Forward used exp(-j*2π*q*carry/L),
            //    so the adjoint here uses the conjugate factor, i.e. exp(+j*2π*q*carry/L).
            var Xhat = new Complex32[Mloc, L];
            var bufL = new Complex32[L];
            var bufL2 = new Complex32[L];

            for (int n0 = 0; n0 < Mloc; n0++)
            {
                // prodMain[q] = FFT_L{ Cmain[n0,:] }
                for (int l = 0; l < L; l++) bufL[l] = Cmain[n0, l];
                FFT(bufL, false);

                // prodHalf[q] = FFT_L{ Chalf[n0,:] }
                for (int l = 0; l < L; l++) bufL2[l] = Chalf[n0, l];
                FFT(bufL2, false);

                int carry = ((n0 + Mloc / 2) >= Mloc) ? 1 : 0;

                for (int q = 0; q < L; q++)
                {
                    // main branch adjoint: multiply by S_hat (adjoint of conj(S_hat))
                    Xhat[n0, q] += S_hat[n0, q] * bufL[q];

                    // half branch adjoint: multiply by T_hat * exp(+j*2π*q*carry/L)
                    if (carry != 0)
                    {
                        float ang = 2f * Maths.Pi * q * carry / L;
                        var shiftPhase = Maths.Exp(Complex32.I * ang); // adjoint to the forward’s (-) sign
                        Xhat[n0, q] += T_hat[n0, q] * shiftPhase * bufL2[q];
                    }
                    else
                    {
                        Xhat[n0, q] += T_hat[n0, q] * bufL2[q];
                    }
                }
            }

            // 3) Return from Zak/FFT domain back to the time domain:
            //    For each residue n0, take IFFT_L over q to get A[r*M + n0].
            var Arec = new Complex32[N];

            for (int n0 = 0; n0 < Mloc; n0++)
            {
                for (int q = 0; q < L; q++) 
                    bufL[q] = Xhat[n0, q];

                FFT(bufL, true); // inverse FFT along r (assumed to divide by L internally)

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
                if (!Maths.IsEven(Mloc)) throw new Exception("M must be even");
                int L = N / Mloc;
                if (L * Mloc != N) throw new Exception("N must be divisible by M");

                // Step 1: Generate the initial WH analysis window g0 of length N
                //         from the current IWindow object (e.g., Gaussian, Hanning, etc.)
                var g0 = WeylHeisenbergTransform.GetPacket(window, N);

                // Step 2: Orthogonalize g0 using Zak-domain orthogonalization
                //         to produce a WH-orthonormal window g.
                //         This matches the behavior of Matrix(..., true) in the slow implementation.
                var zakOrth = new ZakTransform(Mloc);
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
        /// Returns the complex phase factor e^{-j * π * k / 2}.
        /// This term has a period of 4 over k, meaning it only takes
        /// four distinct values: 1, -j, -1, +j.
        /// 
        /// Such factors often appear in WH transforms, DFT twiddle factors,
        /// and modulation/demodulation stages where a quarter-period
        /// phase shift in the complex plane is required.
        /// 
        /// In the Weyl–Heisenberg transform context, this phase multiplier
        /// accounts for the shift by M/4 in the time-frequency tiling,
        /// separating cosine and sine branches during synthesis/analysis.
        /// </summary>
        /// <param name="k">Frequency index (integer).</param>
        /// <returns>Complex value of e^{-j * π * k / 2}.</returns>
        internal static Complex32 PhaseMinusPiOver2(int k)
        {
            // e^{-jπk/2} is periodic with period 4 in k:
            //   k mod 4 = 0 →  1  (0° phase)
            //   k mod 4 = 1 → -j  (-90° phase)
            //   k mod 4 = 2 → -1  (-180° phase)
            //   k mod 4 = 3 → +j  (+90° phase)
            //
            // Bitwise AND with 3 (0b11) is a fast way to compute k % 4.
            return (k & 3) switch
            {
                0 => new Complex32(+1, 0),  // 0°   :  1
                1 => new Complex32(0, -1),  // -90° : -j
                2 => new Complex32(-1, 0),  // 180° : -1
                _ => new Complex32(0, +1),  // +90° : +j
            };
        }

        /// <summary>
        /// Fast Fourier transform.
        /// </summary>
        /// <param name="a">Input</param>
        /// <param name="inverse">Iverse or not</param>
        internal static void FFT(Complex32[] a, bool inverse)
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
