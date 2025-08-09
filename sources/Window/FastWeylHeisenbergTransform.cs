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
    /// Weyl-Heisenberg transforms. This implementation was designed and developed by Valery Asiryan.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastWeylHeisenbergTransform : IWindowTransform, ITransform
    {
        #region Private data
        private IWindow window;
        private int m;
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes fast Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N]</param>
        /// <param name="direction">Processing direction</param>
        public FastWeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical)
        {
            Window = window; M = m; Direction = direction;
        }
        /// <summary>
        /// Gets or sets number of frequency shifts [4, N].
        /// <remarks>
        /// Even number.
        /// </remarks>
        /// </summary>
        public int M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value <= 2 || !Maths.IsEven(value))
                    throw new Exception("Invalid argument value");

                this.m = value;
            }
        }
        /// <summary>
        /// Gets or sets the processing direction.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Gets or sets the window function.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            int Mloc = this.m;
            if (!Maths.IsEven(Mloc)) throw new Exception("M must be even");
            int L = N / Mloc;
            if (L * Mloc != N) throw new Exception("N must be divisible by M");

            // Precompute Zak/FFT caches of the (orthogonalized) window:
            // S_hat[n0, q] = FFT_L { g[r*M + n0] } over r
            // T_hat[n0, q] = FFT_L { g[r*M + (n0+M/2)%M] } over r
            // Dimensions: [M, L]
            lock (locker)
            {
                PrepareFastCachesPolyphase(N, Mloc); // [M, L]
            }

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
            var B = new Complex32[N];
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
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            int Mloc = this.m;
            if (!Maths.IsEven(Mloc)) throw new Exception("M must be even");
            int L = N / Mloc;
            if (L * Mloc != N) throw new Exception("N must be divisible by M");

            // Use the same window caches as in Forward():
            // S_hat/T_hat are FFT_L of the window’s polyphase components.
            lock (locker)
            {
                PrepareFastCachesPolyphase(N, Mloc); // [M, L]
            }

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
                    Cmain[n0, l] = Sp_main[n0] * (1f / Maths.Sqrt(L * N)); // keep your current scale convention
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
                for (int q = 0; q < L; q++) bufL[q] = Xhat[n0, q];
                FFT(bufL, true); // inverse FFT along r (assumed to divide by L internally)

                for (int r = 0; r < L; r++)
                    Arec[r * Mloc + n0] = bufL[r];
            }

            return Arec;
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            Complex32[,] B = (Complex32[,])A.Clone();
            int N = B.GetLength(0), M = B.GetLength(1);

            FastZakTransform zakTransform = new FastZakTransform(m);
            float[] g0 = zakTransform.Orthogonalize(WeylHeisenbergTransform.GetPacket(this.window, N));
            float[] g1 = zakTransform.Orthogonalize(WeylHeisenbergTransform.GetPacket(this.window, M));

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

                    row = Forward(row);

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

                    col = Forward(col);

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

                    col = Forward(col);

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

                    row = Forward(row);

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
        public Complex32[,] Backward(Complex32[,] B)
        {
            Complex32[,] A = (Complex32[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            FastZakTransform zakTransform = new FastZakTransform(m);
            float[] g0 = zakTransform.Orthogonalize(WeylHeisenbergTransform.GetPacket(this.window, N));
            float[] g1 = zakTransform.Orthogonalize(WeylHeisenbergTransform.GetPacket(this.window, M));

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

                    col = Backward(col);

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

                    row = Backward(row);

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

                    col = Backward(col);

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

                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Private voids

        private static readonly FastFourierTransform fastFourierTransform = new FastFourierTransform(false, Direction.Vertical);
        private readonly object locker = new object();
        private Complex32[,] S_hat; // [M, L] FFT_L по r от s_n0[r] = g[r*M + n0]
        private Complex32[,] T_hat; // [M, L] FFT_L по r от t_n0[r] = g[r*M + (n0 + M/2) mod M]
        private int cachedN_poly = -1;

        /// <summary>
        /// Precomputes FFT_L of the polyphase components of the orthogonalized window.
        /// These cached arrays S_hat and T_hat are used in both Forward() and Backward()
        /// to avoid recomputing window FFTs for every transform call.
        /// </summary>
        /// <param name="N">Total signal length.</param>
        /// <param name="Mloc">Number of frequency shifts M (must be even).</param>
        private void PrepareFastCachesPolyphase(int N, int Mloc)
        {
            // If cache for the given N is already computed, reuse it
            if (cachedN_poly == N && S_hat != null && T_hat != null)
                return;

            int L = N / Mloc; // number of time shifts

            // Step 1: Generate the initial WH analysis window g0 of length N
            //         from the current IWindow object (e.g., Gaussian, Hanning, etc.)
            var g0 = WeylHeisenbergTransform.GetPacket(this.window, N);

            // Step 2: Orthogonalize g0 using Zak-domain orthogonalization
            //         to produce a WH-orthonormal window g.
            //         This matches the behavior of Matrix(..., true) in the slow implementation.
            var zakOrth = new ZakTransform(Mloc);
            var g = zakOrth.Orthogonalize(g0); // real-valued array of length N

            // Allocate caches for the FFTs of polyphase components:
            // S_hat[n0, q] = FFT_L over r of polyphase component g[r*M + n0]
            // T_hat[n0, q] = FFT_L over r of polyphase component g[r*M + (n0 + M/2) % M]
            S_hat = new Complex32[Mloc, L];
            T_hat = new Complex32[Mloc, L];

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

            // Mark cache as valid for this N
            cachedN_poly = N;
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
        private static Complex32 PhaseMinusPiOver2(int k)
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

        #endregion
    }
}
