using System;
using UMapx.Core;
using UMapx.Window;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the discrete wavelet.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wavelet
    /// </remarks>
    /// </summary>
    public partial class WaveletPack
    {
        #region Gabor-Zak wavelets
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=4).
        /// </summary>
        public static WaveletPack GZ2
        {
            get
            {
                return WaveletPack.GaborZak(2);
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=8).
        /// </summary>
        public static WaveletPack GZ3
        {
            get
            {
                return WaveletPack.GaborZak(3);
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=16).
        /// </summary>
        public static WaveletPack GZ4
        {
            get
            {
                return WaveletPack.GaborZak(4);
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=32).
        /// </summary>
        public static WaveletPack GZ5
        {
            get
            {
                return WaveletPack.GaborZak(5);
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=64).
        /// </summary>
        public static WaveletPack GZ6
        {
            get
            {
                return WaveletPack.GaborZak(6);
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=128).
        /// </summary>
        public static WaveletPack GZ7
        {
            get
            {
                return WaveletPack.GaborZak(7);
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT).
        /// Total filter length is N = 2^n; the pack contains four N-tap filters
        /// (analysis low/high and synthesis low/high).
        /// </summary>
        /// <param name="n">
        /// Size parameter (2..7): N = 2^n.
        /// </param>
        /// <returns>
        /// WaveletPack with (h0, h1, g0, g1), each of length N = 2^n.
        /// </returns>
        internal static WaveletPack GaborZak(int n)
        {
            if (n < 2) throw new ArgumentException("Wavelet size must be greater or equal 2");
            if (n > 7) throw new ArgumentException("Wavelet size must be less or equal 7");

            // Gabor–Zak wavelets (dyadic, periodic).
            // Builds an orthonormal 2-channel (M=2) wavelet filter bank on ℓ²(ℤ_N)
            // by Zak-domain orthogonalization of a centered truncated Gaussian (Gabor) prototype
            // and CQF construction for the high-pass.
            //
            // Idea: take a symmetric N-periodic Gaussian g → enforce Walnut/Wexler–Raz locally
            // in the Zak domain (|Z₀|² + |Z₁|² = const) → inverse mapping gives the paraunitary
            // low-pass h₀; define h₁ by CQF: h₁[n] = (−1)ⁿ · h₀[Nf−1−n].
            //
            // Properties:
            // • Orthonormal dyadic filter bank (perfect reconstruction, periodic DWT).
            // • Symmetric (N-periodic even) when N ≡ 0 (mod 4) — i.e., n ≥ 2.
            // • Low-pass normalization: Σ h₀ = √2 and H₀(π) = 0 ⇒ ψ has zero mean (≥1 vanishing moment).
            //
            // Notes:
            // • Discrete/periodic construction on ℤ_N (no claim of a continuous MRA on L²(ℝ)).
            // • With N = 2^n and n ≥ 2 we have R = gcd(M, N/2) = 2 (stable two-coset Zak orthogonalization).
            // • This discrete wavelet filter bank was found and introduced by Valery Asiryan (Yerevan, Armenia, 2025).

            var N = (int)Maths.Pow(2, n);
            var window = Gabor.Scaled(frameSize: N);
            var zak = new FastZakTransform(m: 2);
            var g0 = window.GetWindow();
            var h0 = zak.Orthogonalize(g0);
            return WaveletPack.Create(h0);
        }

        #endregion
    }
}
