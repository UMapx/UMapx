using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines a discrete wavelet decomposition.
    /// </summary>
    /// <remarks>
    /// For the correct wavelet transform of a signal, it is necessary that its dimension be a power of 2.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_wavelet_transform
    /// </remarks>
    [Serializable]
    public class WaveletDecomposition : IPyramidTransform
    {
        #region Private data
        private float[] lp, hp, ilp, ihp;
        private bool normalized;
        private int levels;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes a discrete wavelet decomposition.
        /// </summary>
        /// <param name="wavelet">Discrete wavelet</param>
        /// <param name="levels">Number of levels</param>
        /// <param name="normalized">Normalized transform or not</param>
        public WaveletDecomposition(WaveletPacket wavelet, int levels = 1, bool normalized = true)
        {
            Wavelet = wavelet;
            Levels = levels;
            Normalized = normalized;
        }
        /// <summary>
        /// Gets or sets the number of transform levels.
        /// </summary>
        public int Levels
        {
            get
            {
                return levels;
            }
            set
            {
                if (value < 1)
                    throw new ArgumentException("Number of levels cannot be less than 1");

                levels = value;
            }
        }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public bool Normalized
        {
            get 
            {
                return normalized;
            }
            set 
            { 
                normalized = value; 
            }
        }
        /// <summary>
        /// Gets or sets the discrete wavelet.
        /// </summary>
        public WaveletPacket Wavelet
        {
            get
            {
                return new WaveletPacket(lp, hp, ilp, ihp);
            }
            set
            {
                lp = value.LowPass; hp = value.HighPass; ilp = value.ILowPass; ihp = value.IHighPass;
            }
        }

        #endregion

        #region Wavelet decomposition
        /// <summary>
        /// Forward wavelet decomposition.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[][] Forward(float[] A)
        {
            if (A == null) throw new ArgumentNullException(nameof(A));
            int N = A.Length;
            int L = Math.Min((int)Maths.Log2(N), Levels);
            if (L < 1) throw new ArgumentException("Signal length must allow at least 1 level");

            var details = new float[L][];
            float[] a = (float[])A.Clone();

            for (int lev = 1; lev <= L; lev++)
            {
                int bound = a.Length;
                if (!Maths.IsEven(bound)) throw new ArgumentException("Size must be even per level");
                var work = new float[bound];
                DWT1D(a, bound, work);
                int half = bound >> 1;
                var aNext = new float[half];
                var dLev = new float[half];
                Array.Copy(work, 0, aNext, 0, half);
                Array.Copy(work, half, dLev, 0, half);
                details[L - lev] = dLev;
                a = aNext;
            }

            var outArr = new float[L + 1][];
            outArr[0] = a;
            for (int i = 0; i < L; i++) outArr[i + 1] = details[i];
            return outArr;
        }
        /// <summary>
        /// Backward wavelet decomposition.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[][] B)
        {
            if (B == null || B.Length < 2) throw new ArgumentException("Expect at least {A_L, D_L}");
            int L = B.Length - 1;
            float[] a = B[0];

            for (int i = 1; i <= L; i++)
            {
                var d = B[i];
                int half = a.Length; 
                int bound = half << 1; 
                if (!Maths.IsEven(bound)) throw new ArgumentException("Invalid pair lengths");
                // Compose [a|d]
                var a_d = new float[bound];
                Array.Copy(a, 0, a_d, 0, half);
                Array.Copy(d, 0, a_d, half, half);
                var dest = new float[bound];
                IDWT1D(a_d, bound, dest);
                a = dest;
            }
            return a;
        }
        /// <summary>
        /// Forward wavelet decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[][,] Forward(float[,] A)
        {
            if (A == null) throw new ArgumentNullException(nameof(A));
            int R = A.GetLength(0), C = A.GetLength(1);
            int L = Math.Min(Math.Min((int)Maths.Log2(R), (int)Maths.Log2(C)), Levels);
            if (L < 1) throw new ArgumentException("Image size must allow at least 1 level");

            var W = (float[,])A.Clone();
            var packs = new float[1 + 3 * L][,];

            // local working buffers sized to full extent once per call
            var rowBuf = new float[C];
            var colBuf = new float[R];
            var work = new float[Math.Max(R, C)];

            int h = R, w = C;
            for (int lev = 1; lev <= L; lev++)
            {
                int halfW = w >> 1, halfH = h >> 1;
                // rows
                for (int r = 0; r < h; r++)
                {
                    for (int c = 0; c < w; c++) rowBuf[c] = W[r, c];
                    DWT1D(rowBuf, w, work);
                    for (int c = 0; c < halfW; c++) W[r, c] = work[c];
                    for (int c = 0; c < halfW; c++) W[r, halfW + c] = work[halfW + c];
                }
                // columns
                for (int c = 0; c < w; c++)
                {
                    for (int r = 0; r < h; r++) colBuf[r] = W[r, c];
                    DWT1D(colBuf, h, work);
                    for (int r = 0; r < halfH; r++) W[r, c] = work[r];
                    for (int r = 0; r < halfH; r++) W[halfH + r, c] = work[halfH + r];
                }
                // slice details
                var LH = new float[halfH, halfW];
                var HL = new float[halfH, halfW];
                var HH = new float[halfH, halfW];

                for (int r = 0; r < halfH; r++)
                {
                    for (int c = 0; c < halfW; c++)
                    {
                        LH[r, c] = W[r, halfW + c];
                        HL[r, c] = W[halfH + r, c];
                        HH[r, c] = W[halfH + r, halfW + c];
                    }
                }

                packs[1 + 3 * (L - lev) + 0] = LH;
                packs[1 + 3 * (L - lev) + 1] = HL;
                packs[1 + 3 * (L - lev) + 2] = HH;

                h >>= 1; w >>= 1;
            }

            // coarsest LL
            int hL2 = R >> L, wL2 = C >> L;
            var LL = new float[hL2, wL2];
            for (int r = 0; r < hL2; r++) for (int c = 0; c < wL2; c++) LL[r, c] = W[r, c];
            packs[0] = LL;
            return packs;
        }
        /// <summary>
        /// Backward wavelet decomposition.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[][,] B)
        {
            if (B == null || B.Length < 4) throw new ArgumentException("Expect at least {LL_1, LH_1, HL_1, HH_1}");
            if ((B.Length - 1) % 3 != 0) throw new ArgumentException("Invalid pyramid length");
            int L = (B.Length - 1) / 3;
            int h = B[0].GetLength(0), w = B[0].GetLength(1);
            int R = h << L, C = w << L;
            var W = new float[R, C];

            // local working buffers sized to full extent once per call
            var rowBuf = new float[C];
            var colBuf = new float[R];
            var work = new float[Math.Max(R, C)];

            // place coarsest LL
            for (int r = 0; r < h; r++)
            {
                for (int c = 0; c < w; c++)
                {
                    W[r, c] = B[0][r, c];
                }
            }

            for (int lev = L; lev >= 1; lev--)
            {
                int halfH = h, halfW = w; h <<= 1; w <<= 1;
                var LH = B[1 + 3 * (L - lev) + 0];
                var HL = B[1 + 3 * (L - lev) + 1];
                var HH = B[1 + 3 * (L - lev) + 2];

                for (int r = 0; r < halfH; r++)
                {
                    for (int c = 0; c < halfW; c++)
                    {
                        W[r, halfW + c] = LH[r, c];
                        W[halfH + r, c] = HL[r, c];
                        W[halfH + r, halfW + c] = HH[r, c];
                    }
                }

                // columns inverse
                for (int c = 0; c < w; c++)
                {
                    for (int r = 0; r < h; r++) colBuf[r] = W[r, c];
                    Array.Copy(colBuf, 0, work, 0, halfH);
                    Array.Copy(colBuf, halfH, work, halfH, halfH);
                    IDWT1D(work, h, colBuf);
                    for (int r = 0; r < h; r++) W[r, c] = colBuf[r];
                }

                // rows inverse
                for (int r = 0; r < h; r++)
                {
                    for (int c = 0; c < w; c++) rowBuf[c] = W[r, c];
                    Array.Copy(rowBuf, 0, work, 0, halfW);
                    Array.Copy(rowBuf, halfW, work, halfW, halfW);
                    IDWT1D(work, w, rowBuf);
                    for (int c = 0; c < w; c++) W[r, c] = rowBuf[c];
                }
            }
            return W;
        }
        /// <summary>
        /// Forward wavelet decomposition.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[][] Forward(Complex32[] A)
        {
            if (A == null) throw new ArgumentNullException(nameof(A));
            int N = A.Length;
            int L = Math.Min((int)Maths.Log2(N), Levels);
            if (L < 1) throw new ArgumentException("Signal length must allow at least 1 level");
            var details = new Complex32[L][];
            Complex32[] a = (Complex32[])A.Clone();

            for (int lev = 1; lev <= L; lev++)
            {
                int bound = a.Length; if (!Maths.IsEven(bound)) throw new ArgumentException("Size must be even per level");
                var work = new Complex32[bound];
                DWT1D(a, bound, work);
                int half = bound >> 1;
                var aNext = new Complex32[half];
                var dLev = new Complex32[half];
                Array.Copy(work, 0, aNext, 0, half);
                Array.Copy(work, half, dLev, 0, half);
                details[L - lev] = dLev;
                a = aNext;
            }

            var outArr = new Complex32[L + 1][];
            outArr[0] = a;
            for (int i = 0; i < L; i++) outArr[i + 1] = details[i];
            return outArr;
        }
        /// <summary>
        /// Backward wavelet decomposition.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[][] B)
        {
            if (B == null || B.Length < 2) throw new ArgumentException("Expect at least {A_L, D_L}");
            int L = B.Length - 1;
            Complex32[] a = B[0];

            for (int i = 1; i <= L; i++)
            {
                var d = B[i];
                int half = a.Length; 
                int bound = half << 1; 
                if (!Maths.IsEven(bound)) throw new ArgumentException("Invalid pair lengths");
                // Compose [a|d]
                var a_d = new Complex32[bound];
                Array.Copy(a, 0, a_d, 0, half);
                Array.Copy(d, 0, a_d, half, half);
                var dest = new Complex32[bound];
                IDWT1D(a_d, bound, dest);
                a = dest;
            }
            return a;
        }
        /// <summary>
        /// Forward wavelet decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[][,] Forward(Complex32[,] A)
        {
            if (A == null) throw new ArgumentNullException(nameof(A));
            int R = A.GetLength(0), C = A.GetLength(1);
            int L = Math.Min(Math.Min((int)Maths.Log2(R), (int)Maths.Log2(C)), Levels);
            if (L < 1) throw new ArgumentException("Image size must allow at least 1 level");

            var W = (Complex32[,])A.Clone();
            var packs = new Complex32[1 + 3 * L][,];

            // local working buffers sized to full extent once per call
            var rowBuf = new Complex32[C];
            var colBuf = new Complex32[R];
            var work = new Complex32[Math.Max(R, C)];

            int h = R, w = C;
            for (int lev = 1; lev <= L; lev++)
            {
                int halfW = w >> 1, halfH = h >> 1;
                // rows
                for (int r = 0; r < h; r++)
                {
                    for (int c = 0; c < w; c++) rowBuf[c] = W[r, c];
                    DWT1D(rowBuf, w, work);
                    for (int c = 0; c < halfW; c++) W[r, c] = work[c];
                    for (int c = 0; c < halfW; c++) W[r, halfW + c] = work[halfW + c];
                }
                // columns
                for (int c = 0; c < w; c++)
                {
                    for (int r = 0; r < h; r++) colBuf[r] = W[r, c];
                    DWT1D(colBuf, h, work);
                    for (int r = 0; r < halfH; r++) W[r, c] = work[r];
                    for (int r = 0; r < halfH; r++) W[halfH + r, c] = work[halfH + r];
                }
                // slice details
                var LH = new Complex32[halfH, halfW];
                var HL = new Complex32[halfH, halfW];
                var HH = new Complex32[halfH, halfW];
                for (int r = 0; r < halfH; r++)
                {
                    for (int c = 0; c < halfW; c++)
                    { 
                        LH[r, c] = W[r, halfW + c]; 
                        HL[r, c] = W[halfH + r, c];
                        HH[r, c] = W[halfH + r, halfW + c]; 
                    }
                }
                packs[1 + 3 * (L - lev) + 0] = LH;
                packs[1 + 3 * (L - lev) + 1] = HL;
                packs[1 + 3 * (L - lev) + 2] = HH;

                h >>= 1; w >>= 1;
            }

            // coarsest LL
            int hL2 = R >> L, wL2 = C >> L;
            var LL = new Complex32[hL2, wL2];

            for (int r = 0; r < hL2; r++)
            {
                for (int c = 0; c < wL2; c++)
                {
                    LL[r, c] = W[r, c];
                }
            }
            packs[0] = LL;
            return packs;
        }
        /// <summary>
        /// Backward wavelet decomposition.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[][,] B)
        {
            if (B == null || B.Length < 4) throw new ArgumentException("Expect at least {LL_1, LH_1, HL_1, HH_1}");
            if ((B.Length - 1) % 3 != 0) throw new ArgumentException("Invalid pyramid length");
            int L = (B.Length - 1) / 3;
            int h = B[0].GetLength(0), w = B[0].GetLength(1);
            int R = h << L, C = w << L;
            var W = new Complex32[R, C];

            // local working buffers sized to full extent once per call
            var rowBuf = new Complex32[C];
            var colBuf = new Complex32[R];
            var work = new Complex32[Math.Max(R, C)];

            // place coarsest LL
            for (int r = 0; r < h; r++)
            {
                for (int c = 0; c < w; c++)
                {
                    W[r, c] = B[0][r, c];
                }
            }

            for (int lev = L; lev >= 1; lev--)
            {
                int halfH = h, halfW = w; h <<= 1; w <<= 1;
                var LH = B[1 + 3 * (L - lev) + 0];
                var HL = B[1 + 3 * (L - lev) + 1];
                var HH = B[1 + 3 * (L - lev) + 2];

                for (int r = 0; r < halfH; r++)
                {
                    for (int c = 0; c < halfW; c++)
                    {
                        W[r, halfW + c] = LH[r, c];
                        W[halfH + r, c] = HL[r, c];
                        W[halfH + r, halfW + c] = HH[r, c];
                    }
                }

                // columns inverse
                for (int c = 0; c < w; c++)
                {
                    for (int r = 0; r < h; r++) colBuf[r] = W[r, c];
                    Array.Copy(colBuf, 0, work, 0, halfH);
                    Array.Copy(colBuf, halfH, work, halfH, halfH);
                    IDWT1D(work, h, colBuf);
                    for (int r = 0; r < h; r++) W[r, c] = colBuf[r];
                }

                // rows inverse
                for (int r = 0; r < h; r++)
                {
                    for (int c = 0; c < w; c++) rowBuf[c] = W[r, c];
                    Array.Copy(rowBuf, 0, work, 0, halfW);
                    Array.Copy(rowBuf, halfW, work, halfW, halfW);
                    IDWT1D(work, w, rowBuf);
                    for (int c = 0; c < w; c++) W[r, c] = rowBuf[c];
                }
            }
            return W;
        }
        #endregion

        #region Private methods
        /// <summary>
        /// Performs a single-level 1D wavelet <b>analysis</b> (decimation by 2) with periodic boundaries.
        /// Writes concatenated coefficients into <paramref name="output"/> as <c>[A(0..h-1), D(0..h-1)]</c>,
        /// where <c>h = bound/2</c>.
        /// </summary>
        /// <param name="input">Source samples; only the first <paramref name="bound"/> values are used</param>
        /// <param name="bound">Working length (must be even). Defines <c>h = bound/2</c></param>
        /// <param name="output">Destination buffer of length <paramref name="bound"/> receiving A then D</param>
        /// <remarks>
        /// Uses circular (periodic) extension by advancing rotating indices instead of using modulo per tap.
        /// If <see cref="Normalized"/> is true, the output bands are scaled by <c>1/√2</c> to match orthonormal energy.
        /// Complexity: O(bound · (|lp| + |hp|)).
        /// </remarks>
        private void DWT1D(float[] input, int bound, float[] output)
        {
            if (!Maths.IsEven(bound)) bound--;
            int h = bound >> 1;
            int lpLen = lp.Length, hpLen = hp.Length;
            int lpStart = -((lpLen >> 1) - 1);
            int hpStart = -((hpLen >> 1) - 1);
            int cL0 = ModBound(lpStart, bound);
            int cH0 = ModBound(hpStart, bound);

            for (int r = 0; r < h; r++)
            {
                float a = 0f, b = 0f; int cL = cL0, cH = cH0;
                for (int k = 0; k < lpLen; k++) { a += lp[k] * input[cL]; if (++cL == bound) cL = 0; }
                for (int k = 0; k < hpLen; k++) { b += hp[k] * input[cH]; if (++cH == bound) cH = 0; }
                cL0 += 2; if (cL0 >= bound) cL0 -= bound; cH0 += 2; if (cH0 >= bound) cH0 -= bound;
                if (normalized) { output[r] = a / Maths.Sqrt2; output[r + h] = b / Maths.Sqrt2; }
                else { output[r] = a; output[r + h] = b; }
            }
        }
        /// <summary>
        /// Performs a single-level 1D wavelet <b>synthesis</b> from concatenated bands.
        /// Reconstructs a signal of length <paramref name="bound"/> from <c>a_d = [A(0..h-1), D(0..h-1)]</c>, <c>h = bound/2</c>.
        /// </summary>
        /// <param name="a_d">Input buffer holding A followed by D coefficients</param>
        /// <param name="bound">Output length (must be even). Defines <c>h = bound/2</c></param>
        /// <param name="dest">Destination signal of length <paramref name="bound"/></param>
        /// <remarks>
        /// Uses odd-phase upsampling (values placed at indices <c>i+1</c>) for both A and D branches to reproduce legacy phasing.
        /// Circular (periodic) extension is applied during convolution; indices advance by one per output sample.
        /// If <see cref="Normalized"/> is true, the result is scaled by <c>√2</c> (inverse of analysis scaling).
        /// Complexity: O(bound · (|ilp| + |ihp|)).
        /// </remarks>
        private void IDWT1D(float[] a_d, int bound, float[] dest)
        {
            if (!Maths.IsEven(bound)) bound--;
            int h = bound >> 1;
            var upL = new float[bound];
            var upH = new float[bound];

            Array.Clear(upL, 0, bound); Array.Clear(upH, 0, bound);
            for (int j = 0, i = 0; j < h; j++, i += 2) { upL[i + 1] = a_d[j]; upH[i + 1] = a_d[h + j]; }

            int lpLen = ilp.Length, hpLen = ihp.Length;
            int lpStart = -((lpLen >> 1) - 1);
            int hpStart = -((hpLen >> 1) - 1);
            int cL0 = ModBound(lpStart, bound);
            int cH0 = ModBound(hpStart, bound);

            for (int i = 0; i < bound; i++)
            {
                float s = 0f; int cL = cL0, cH = cH0;
                for (int k = 0; k < lpLen; k++) { s += ilp[k] * upL[cL]; if (++cL == bound) cL = 0; }
                for (int k = 0; k < hpLen; k++) { s += ihp[k] * upH[cH]; if (++cH == bound) cH = 0; }
                dest[i] = normalized ? s * Maths.Sqrt2 : s;
                cL0++; if (cL0 == bound) cL0 = 0; cH0++; if (cH0 == bound) cH0 = 0;
            }
        }
        /// <summary>
        /// Performs a single-level 1D wavelet <b>analysis</b> (decimation by 2) with periodic boundaries.
        /// Writes concatenated coefficients into <paramref name="output"/> as <c>[A(0..h-1), D(0..h-1)]</c>,
        /// where <c>h = bound/2</c>.
        /// </summary>
        /// <param name="input">Source samples; only the first <paramref name="bound"/> values are used</param>
        /// <param name="bound">Working length (must be even). Defines <c>h = bound/2</c></param>
        /// <param name="output">Destination buffer of length <paramref name="bound"/> receiving A then D</param>
        /// <remarks>
        /// Uses circular (periodic) extension by advancing rotating indices instead of using modulo per tap.
        /// If <see cref="Normalized"/> is true, the output bands are scaled by <c>1/√2</c> to match orthonormal energy.
        /// Complexity: O(bound · (|lp| + |hp|)).
        /// </remarks>
        private void DWT1D(Complex32[] input, int bound, Complex32[] output)
        {
            if (!Maths.IsEven(bound)) bound--;
            int h = bound >> 1;
            int lpLen = lp.Length, hpLen = hp.Length;
            int lpStart = -((lpLen >> 1) - 1);
            int hpStart = -((hpLen >> 1) - 1);
            int cL0 = ModBound(lpStart, bound);
            int cH0 = ModBound(hpStart, bound);
            for (int r = 0; r < h; r++)
            {
                Complex32 a = 0, b = 0; int cL = cL0, cH = cH0;
                for (int k = 0; k < lpLen; k++) { a += lp[k] * input[cL]; if (++cL == bound) cL = 0; }
                for (int k = 0; k < hpLen; k++) { b += hp[k] * input[cH]; if (++cH == bound) cH = 0; }
                cL0 += 2; if (cL0 >= bound) cL0 -= bound; cH0 += 2; if (cH0 >= bound) cH0 -= bound;
                output[r] = normalized ? a / Maths.Sqrt2 : a;
                output[r + h] = normalized ? b / Maths.Sqrt2 : b;
            }
        }
        /// <summary>
        /// Performs a single-level 1D wavelet <b>synthesis</b> from concatenated bands.
        /// Reconstructs a signal of length <paramref name="bound"/> from <c>a_d = [A(0..h-1), D(0..h-1)]</c>, <c>h = bound/2</c>.
        /// </summary>
        /// <param name="a_d">Input buffer holding A followed by D coefficients</param>
        /// <param name="bound">Output length (must be even). Defines <c>h = bound/2</c></param>
        /// <param name="dest">Destination signal of length <paramref name="bound"/></param>
        /// <remarks>
        /// Uses odd-phase upsampling (values placed at indices <c>i+1</c>) for both A and D branches to reproduce legacy phasing.
        /// Circular (periodic) extension is applied during convolution; indices advance by one per output sample.
        /// If <see cref="Normalized"/> is true, the result is scaled by <c>√2</c> (inverse of analysis scaling).
        /// Complexity: O(bound · (|ilp| + |ihp|)).
        /// </remarks>
        private void IDWT1D(Complex32[] a_d, int bound, Complex32[] dest)
        {
            if (!Maths.IsEven(bound)) bound--;
            int h = bound >> 1;
            var upL = new Complex32[bound];
            var upH = new Complex32[bound];

            Array.Clear(upL, 0, bound); Array.Clear(upH, 0, bound);
            for (int j = 0, i = 0; j < h; j++, i += 2) { upL[i + 1] = a_d[j]; upH[i + 1] = a_d[h + j]; }
            int lpLen = ilp.Length, hpLen = ihp.Length;
            int lpStart = -((lpLen >> 1) - 1); int hpStart = -((hpLen >> 1) - 1);
            int cL0 = ModBound(lpStart, bound); int cH0 = ModBound(hpStart, bound);
            for (int i = 0; i < bound; i++)
            {
                Complex32 s = 0; int cL = cL0, cH = cH0;
                for (int k = 0; k < lpLen; k++) { s += ilp[k] * upL[cL]; if (++cL == bound) cL = 0; }
                for (int k = 0; k < hpLen; k++) { s += ihp[k] * upH[cH]; if (++cH == bound) cH = 0; }
                dest[i] = normalized ? s * Maths.Sqrt2 : s;
                cL0++; if (cL0 == bound) cL0 = 0; cH0++; if (cH0 == bound) cH0 = 0;
            }
        }
        /// <summary>
        /// Returns <paramref name="j"/> wrapped into the range <c>[0, n)</c>.
        /// </summary>
        /// <param name="j">Index (may be negative or ≥ <paramref name="n"/>)</param>
        /// <param name="n">Modulus (&gt; 0)</param>
        /// <returns>Value in <c>[0, n)</c> equivalent to <paramref name="j"/> modulo <paramref name="n"/></returns>
        private static int ModBound(int j, int n) { int r = j % n; return r < 0 ? r + n : r; }
        #endregion
    }
}