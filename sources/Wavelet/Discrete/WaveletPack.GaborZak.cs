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
                // Original procedure:
                //return WaveletPack.GaborZak(2);
                return WaveletPack.Create(new float[] {
                    0f,
                    0.70710677f,
                    0.70710677f,
                    0f
                });
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=8).
        /// </summary>
        public static WaveletPack GZ3
        {
            get
            {
                // Original procedure:
                //return WaveletPack.GaborZak(3);
                return WaveletPack.Create(new float[] {
                    0.0012465268f,
                    -0.029662542f,
                    0.02966249f,
                    0.7058603f,
                    0.7058603f,
                    0.029662471f,
                    -0.02966255f,
                    0.0012465119f
                });
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=16).
        /// </summary>
        public static WaveletPack GZ4
        {
            get
            {
                // Original procedure:
                //return WaveletPack.GaborZak(4);
                return WaveletPack.Create(new float[] {
                    0.001059331f,
                    -0.0049649402f,
                    0.00089050457f,
                    0.02053153f,
                    -0.008592626f,
                    -0.09199238f,
                    0.09606685f,
                    0.6941085f,
                    0.69410855f,
                    0.09606686f,
                    -0.09199238f,
                    -0.008592622f,
                    0.02053153f,
                    0.0008904934f,
                    -0.0049649477f,
                    0.0010593459f
                });
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=32).
        /// </summary>
        public static WaveletPack GZ5
        {
            get
            {
                // Original procedure:
                //return WaveletPack.GaborZak(5);
                return WaveletPack.Create(new float[] {
                    0.0006670952f,
                    -0.0014538579f,
                    -0.00018056482f,
                    0.0029585995f,
                    -0.00029181875f,
                    -0.0059948377f,
                    0.0011297036f,
                    0.012299581f,
                    -0.0032829107f,
                    -0.025798228f,
                    0.009965725f,
                    0.056083903f,
                    -0.034260236f,
                    -0.13308713f,
                    0.15541902f,
                    0.67293245f,
                    0.6729324f,
                    0.15541911f,
                    -0.13308713f,
                    -0.03426018f,
                    0.05608394f,
                    0.009965759f,
                    -0.025798174f,
                    -0.0032828925f,
                    0.012299588f,
                    0.001129698f,
                    -0.005994808f,
                    -0.00029180385f,
                    0.0029586367f,
                    -0.00018055737f,
                    -0.0014538504f,
                    0.00066710636f
                });
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=64).
        /// </summary>
        public static WaveletPack GZ6
        {
            get
            {
                // Original procedure:
                //return WaveletPack.GaborZak(6);
                return WaveletPack.Create(new float[] {
                    0.0003711935f,
                    -0.0005511902f,
                    -0.00023291633f,
                    0.0007955674f,
                    0.00011730194f,
                    -0.0011357423f,
                    -8.72463E-06f,
                    0.0016153548f,
                    -0.00011170376f,
                    -0.0022981297f,
                    0.00026626047f,
                    0.0032761283f,
                    -0.0004915516f,
                    -0.004687311f,
                    0.0008459219f,
                    0.0067367246f,
                    -0.0014368212f,
                    -0.00973925f,
                    0.0024624513f,
                    0.014180105f,
                    -0.0043130163f,
                    -0.020838853f,
                    0.0077769845f,
                    0.031022813f,
                    -0.014545854f,
                    -0.047203172f,
                    0.028619789f,
                    0.07537001f,
                    -0.062068027f,
                    -0.139644f,
                    0.1863665f,
                    0.65658665f,
                    0.65658706f,
                    0.18636644f,
                    -0.13964348f,
                    -0.062067937f,
                    0.0753705f,
                    0.028620087f,
                    -0.04720272f,
                    -0.014545446f,
                    0.031023262f,
                    0.007777346f,
                    -0.02083858f,
                    -0.0043127704f,
                    0.01418031f,
                    0.0024627056f,
                    -0.009738997f,
                    -0.0014365532f,
                    0.0067368997f,
                    0.00084614847f,
                    -0.0046872394f,
                    -0.00049137f,
                    0.0032761414f,
                    0.0002664579f,
                    -0.0022980953f,
                    -0.0001115324f,
                    0.0016153418f,
                    -8.6613E-06f,
                    -0.0011357255f,
                    0.00011746399f,
                    0.0007954668f,
                    -0.00023279525f,
                    -0.00055134855f,
                    0.00037127547f
                });
            }
        }
        /// <summary>
        /// Returns a symmetric dyadic Gabor–Zak wavelet filter bank (periodic DWT for N=128).
        /// </summary>
        public static WaveletPack GZ7
        {
            get
            {
                // Original procedure:
                //return WaveletPack.GaborZak(7);
                return WaveletPack.Create(new float[] {
                    0.0001938641f,
                    -0.00023879297f,
                    -0.00015845615f,
                    0.00028761942f,
                    0.00012333132f,
                    -0.0003466187f,
                    -9.512529E-05f,
                    0.00041576754f,
                    6.3076615E-05f,
                    -0.0004964685f,
                    -3.932137E-05f,
                    0.0005912734f,
                    1.09523535E-05f,
                    -0.00070687477f,
                    1.4455989E-05f,
                    0.0008402737f,
                    -5.004648E-05f,
                    -0.0010063332f,
                    8.182041E-05f,
                    0.0011977609f,
                    -0.00012414716f,
                    -0.0014316873f,
                    0.00016986812f,
                    0.0017086524f,
                    -0.00023080129f,
                    -0.0020461385f,
                    0.0003019222f,
                    0.0024459509f,
                    -0.00039664295f,
                    -0.002933499f,
                    0.0005151505f,
                    0.0035173732f,
                    -0.0006772806f,
                    -0.0042354036f,
                    0.00087627093f,
                    0.00509063f,
                    -0.0011532516f,
                    -0.006144617f,
                    0.0015166758f,
                    0.007422122f,
                    -0.0020132181f,
                    -0.008998428f,
                    0.0026811808f,
                    0.010927156f,
                    -0.0036102221f,
                    -0.013325617f,
                    0.0048949453f,
                    0.016302515f,
                    -0.0067174304f,
                    -0.020082403f,
                    0.009310857f,
                    0.024914667f,
                    -0.013110357f,
                    -0.031320505f,
                    0.0188319f,
                    0.04019787f,
                    -0.027923234f,
                    -0.05363105f,
                    0.043845765f,
                    0.077370994f,
                    -0.07776834f,
                    -0.13591355f,
                    0.20025751f,
                    0.64701927f,
                    0.6470252f,
                    0.20025848f,
                    -0.13590756f,
                    -0.0777637f,
                    0.077375464f,
                    0.043845847f,
                    -0.053625546f,
                    -0.027920969f,
                    0.04020568f,
                    0.018834446f,
                    -0.031314325f,
                    -0.013103597f,
                    0.024918739f,
                    0.0093167415f,
                    -0.020073134f,
                    -0.0067110816f,
                    0.016311014f,
                    0.004905369f,
                    -0.013320625f,
                    -0.0036070403f,
                    0.010931967f,
                    0.0026845052f,
                    -0.008997273f,
                    -0.0020089583f,
                    0.007426967f,
                    0.0015181f,
                    -0.006143208f,
                    -0.0011501594f,
                    0.0050945487f,
                    0.000879618f,
                    -0.0042304737f,
                    -0.0006701682f,
                    0.0035206988f,
                    0.0005184032f,
                    -0.0029324645f,
                    -0.000396105f,
                    0.0024475143f,
                    0.0003043809f,
                    -0.0020444507f,
                    -0.00022883993f,
                    0.0017088358f,
                    0.00017015124f,
                    -0.0014307657f,
                    -0.00012258813f,
                    0.0011970378f,
                    8.3141495E-05f,
                    -0.0010029036f,
                    -4.6831556E-05f,
                    0.00084154215f,
                    1.857616E-05f,
                    -0.0007071579f,
                    1.0117888E-05f,
                    0.0005931789f,
                    -3.6641024E-05f,
                    -0.000498306f,
                    6.5766275E-05f,
                    0.00041304715f,
                    -9.3222596E-05f,
                    -0.00034669973f,
                    0.00012325775f,
                    0.0002864981f,
                    -0.0001573367f,
                    -0.00023803674f,
                    0.00019518659f
                });
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

            var N = (int)MathF.Pow(2, n);
            var window = Gabor.Scaled(frameSize: N);
            var zak = new FastZakTransform(m: 2);
            var g0 = window.GetWindow();
            var h0 = zak.Orthogonalize(g0);
            return WaveletPack.Create(h0);
        }

        #endregion
    }
}
