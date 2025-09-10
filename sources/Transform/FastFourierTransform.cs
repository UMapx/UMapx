using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Fourier transform using the Cooley-Tukey and Bluestein algorithms.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
    /// https://en.wikipedia.org/wiki/Chirp_Z-transform
    /// </remarks>
    [Serializable]
    public class FastFourierTransform : TransformBaseComplex32, ITransform
    {
        #region Private data
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        private bool normalized;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Fourier transform using the Cooley-Tukey and Bluestein algorithms.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastFourierTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.Normalized = normalized;
            this.Direction = direction;
        }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.normalized;
            }
            set
            {
                this.normalized = value;
            }
        }
        #endregion

        #region Fast Fourier Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            var B = (Complex32[])A.Clone();

            if ((N & (N - 1)) == 0)
                CooleyTukeyFFT(B, inverse: false);
            else
                BluesteinFFT(B, inverse: false);

            if (normalized)
                Scale(B, 1f / Maths.Sqrt(N));

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            var A = (Complex32[])B.Clone();

            if ((N & (N - 1)) == 0)
                CooleyTukeyFFT(A, inverse: true);
            else
                BluesteinFFT(A, inverse: true);

            if (normalized)
                Scale(A, 1f / Maths.Sqrt(N));
            //else
            //    Scale(A, 1f / N);

            return A;
        }
        #endregion

        #region Core FFTs
        /// <summary>
        /// Fast Fourier transform (Bluestein FFT).
        /// </summary>
        /// <param name="data">Array</param>
        /// <param name="inverse">Inverse or not</param>
        private static void BluesteinFFT(Complex32[] data, bool inverse)
        {
            if (data == null || data.Length <= 1) return;

            if (inverse) SwapRealImagInPlace(data);   // conj trick
            TransformBluestein(data);                 // forward kernel
            if (inverse) SwapRealImagInPlace(data);   // back from conj
        }
        /// <summary>
        /// Fast Fourier transform (Cooley-Tukey FFT).
        /// </summary>
        /// <param name="data">Array</param>
        /// <param name="inverse">Inverse or not</param>
        private static void CooleyTukeyFFT(Complex32[] data, bool inverse)
        {
            if (data == null || data.Length <= 1) return;
            int n = data.Length;
            if ((n & (n - 1)) != 0)
                throw new ArgumentException("Length must be a power of 2");

            if (inverse) SwapRealImagInPlace(data);   // conj trick
            TransformRadix2(data);                    // forward kernel
            if (inverse) SwapRealImagInPlace(data);   // back from conj
        }
        /// <summary>
        /// Radix-2 forward kernel (no scaling). Assumes length is a power of two.
        /// </summary>
        /// <param name="data">Array</param>
        /// <exception cref="ArgumentException">Exception</exception>
        private static void TransformRadix2(Complex32[] data)
        {
            int n = data.Length;
            int levels = FloorLog2(n);
            if ((1 << levels) != n) throw new ArgumentException("Length is not a power of 2");

            // trig tables (size n/2): cos(2π k / n), sin(2π k / n)
            var ct = CosTable(n >> 1);
            var st = SinTable(n >> 1);

            // bit-reversed permutation
            for (int i = 0; i < n; i++)
            {
                int j = (int)((uint)ReverseBits32(i) >> (32 - levels));
                if (j > i)
                {
                    var tmp = data[i]; data[i] = data[j]; data[j] = tmp;
                }
            }

            // butterflies
            for (int size = 2; size <= n; size <<= 1)
            {
                int half = size >> 1;
                int step = n / size;
                for (int i = 0; i < n; i += size)
                {
                    int k = 0;
                    for (int j = i; j < i + half; j++, k += step)
                    {
                        int h = j + half;

                        double re = data[h].Real;
                        double im = data[h].Imag;

                        double tpre = +re * ct[k] + im * st[k];
                        double tpim = -re * st[k] + im * ct[k];

                        double rej = data[j].Real;
                        double imj = data[j].Imag;

                        data[h] = new Complex32((float)(rej - tpre), (float)(imj - tpim));
                        data[j] = new Complex32((float)(rej + tpre), (float)(imj + tpim));
                    }
                }
                if (size == n) break; // prevent overflow
            }
        }
        /// <summary>
        /// Bluestein forward kernel (no scaling). Works for arbitrary lengths.
        /// </summary>
        /// <param name="data">Array</param>
        private static void TransformBluestein(Complex32[] data)
        {
            int n = data.Length;
            int m = HighestOneBit(n * 2 + 1) << 1; // next power-of-two >= 2n

            var expC = ExpCosTable(n);
            var expS = ExpSinTable(n);

            // A(z) = x * w1 pre-chirp
            var are = new double[m];
            var aim = new double[m];
            for (int i = 0; i < n; i++)
            {
                double re = data[i].Real, im = data[i].Imag;
                are[i] = +re * expC[i] + im * expS[i];
                aim[i] = -re * expS[i] + im * expC[i];
            }

            // B(z) = chirp kernel (symmetric)
            var bre = new double[m];
            var bim = new double[m];
            bre[0] = expC[0]; bim[0] = expS[0];
            for (int i = 1; i < n; i++)
            {
                bre[i] = bre[m - i] = expC[i];
                bim[i] = bim[m - i] = expS[i];
            }

            // C = A ⊛ B (circular convolution via FFT)
            var cre = new double[m];
            var cim = new double[m];
            Convolve(are, aim, bre, bim, cre, cim);

            // post-chirp: y = C * w1
            for (int i = 0; i < n; i++)
            {
                double re = +cre[i] * expC[i] + cim[i] * expS[i];
                double im = -cre[i] * expS[i] + cim[i] * expC[i];
                data[i] = new Complex32((float)re, (float)im);
            }
        }
        /// <summary>
        /// Circular convolution using in-class FFT kernels. Applies final 1/n scaling.
        /// </summary>
        /// <param name="xre">X.Re</param>
        /// <param name="xim">X.Im</param>
        /// <param name="yre">Y.Re</param>
        /// <param name="yim">Y.Im</param>
        /// <param name="ore">O.Re</param>
        /// <param name="oim">O.Im</param>
        private static void Convolve(double[] xre, double[] xim, double[] yre, double[] yim, double[] ore, double[] oim)
        {
            int n = xre.Length;

            // pack to Complex
            var X = new Complex32[n];
            var Y = new Complex32[n];

            for (int i = 0; i < n; i++)
            {
                X[i] = new Complex32((float)xre[i], (float)xim[i]);
                Y[i] = new Complex32((float)yre[i], (float)yim[i]);
            }

            // FFT forward
            if ((n & (n - 1)) == 0)
            {
                TransformRadix2(X);
                TransformRadix2(Y);
            }
            else
            {
                TransformBluestein(X);
                TransformBluestein(Y);
            }

            // pointwise multiply
            for (int i = 0; i < n; i++)
            {
                float ar = X[i].Real, ai = X[i].Imag;
                float br = Y[i].Real, bi = Y[i].Imag;
                X[i] = new Complex32(ar * br - ai * bi, ai * br + ar * bi);
            }

            // IDFT without scaling (swap trick)
            SwapRealImagInPlace(X);
            if ((n & (n - 1)) == 0) TransformRadix2(X); else TransformBluestein(X);
            SwapRealImagInPlace(X);

            // final scaling by 1/n
            float s = 1f / n;
            for (int i = 0; i < n; i++)
            {
                ore[i] = X[i].Real * s;
                oim[i] = X[i].Imag * s;
            }
        }
        /// <summary>
        /// Swaps real and imaginary parts in-place.
        /// This is equivalent to applying the conjugation trick to switch between forward and inverse.
        /// </summary>
        /// <param name="a">Array</param>
        private static void SwapRealImagInPlace(Complex32[] a)
        {
            for (int i = 0; i < a.Length; i++)
            {
                float r = a[i].Real, im = a[i].Imag;
                a[i] = new Complex32(im, r);
            }
        }
        /// <summary>
        /// Scales a complex vector by a scalar (in-place).
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="s">Scale param</param>
        private static void Scale(Complex32[] a, float s)
        {
            for (int i = 0; i < a.Length; i++)
            {
                a[i].Real *= s; a[i].Imag *= s;
            }
        }
        #endregion

        #region Helpers: tables, bit tricks

        /// <summary>
        /// Returns the highest one bit of <paramref name="i"/> (as a power-of-two integer).
        /// </summary>
        /// <param name="i">Value</param>
        /// <returns>Value</returns>
        private static int HighestOneBit(int i)
        {
            i |= (i >> 1);
            i |= (i >> 2);
            i |= (i >> 4);
            i |= (i >> 8);
            i |= (i >> 16);
            return i - (int)((uint)i >> 1);
        }
        /// <summary>
        /// Reverses the bits of a 32-bit integer.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        private static uint ReverseBits32(int x)
        {
            uint i = (uint)x;
            i = (i & 0x55555555u) << 1 | (i >> 1) & 0x55555555u;
            i = (i & 0x33333333u) << 2 | (i >> 2) & 0x33333333u;
            i = (i & 0x0f0f0f0fu) << 4 | (i >> 4) & 0x0f0f0f0fu;
            i = (i << 24) | ((i & 0xff00u) << 8) | ((i >> 8) & 0xff00u) | (i >> 24);
            return i;
        }
        /// <summary>
        /// Returns floor(log2(n)) for n &gt; 0.
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Value</returns>
        private static int FloorLog2(int n)
        {
            int r = 0; while ((1 << r) < n) r++; return r - (((1 << r) == n) ? 0 : 1);
        }
        /// <summary>
        /// Half-period cosine table: cos(2πk/n) for k=0..n/2-1.
        /// </summary>
        /// <param name="halfN">Value</param>
        /// <returns>Array</returns>
        private static double[] CosTable(int halfN)
        {
            var cosTable = new double[halfN];
            for (int i = 0; i < halfN; i++) cosTable[i] = Math.Cos(Math.PI * i / halfN);
            return cosTable;
        }
        /// <summary>
        /// Half-period sine table: sin(2πk/n) for k=0..n/2-1.
        /// </summary>
        /// <param name="halfN">Value</param>
        /// <returns>Array</returns>
        private static double[] SinTable(int halfN)
        {
            var sinTable = new double[halfN];
            for (int i = 0; i < halfN; i++) sinTable[i] = Math.Sin(Math.PI * i / halfN);
            return sinTable;
        }
        /// <summary>
        /// Chirp cosine table: cos(π i² / n), i=0..n-1 (mod 2n to reduce overflow/rounding).
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Array</returns>
        private static double[] ExpCosTable(int n)
        {
            var expCosTable = new double[n];
            for (int i = 0; i < n; i++)
            {
                int j = (int)((long)i * i % (2L * n));
                expCosTable[i] = Math.Cos(Math.PI * j / n);
            }
            return expCosTable;
        }
        /// <summary>
        /// Chirp sine table: sin(π i² / n), i=0..n-1 (mod 2n to reduce overflow/rounding).
        /// </summary>
        /// <param name="n">Value</param>
        /// <returns>Array</returns>
        private static double[] ExpSinTable(int n)
        {
            var expSinTable = new double[n];
            for (int i = 0; i < n; i++)
            {
                int j = (int)((long)i * i % (2L * n));
                expSinTable[i] = Math.Sin(Math.PI * j / n);
            }
            return expSinTable;
        }
        #endregion
    }
}
