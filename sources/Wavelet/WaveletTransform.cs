using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines a discrete wavelet transform.
    /// </summary>
    /// <remarks>
    /// For the correct wavelet transform of a signal, it is necessary that its dimension be a power of 2.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_wavelet_transform
    /// </remarks>
    [Serializable]
    public class WaveletTransform : IWaveletTransform, ITransform
    {
        #region Private data
        private WaveletDecomposition waveletDecomposition;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the discrete wavelet transform.
        /// </summary>
        /// <param name="waveletDecomposition">Discrete wavelet decomposition</param>
        public WaveletTransform(WaveletDecomposition waveletDecomposition)
        {
            WaveletDecomposition = waveletDecomposition;
        }
        /// <summary>
        /// Gets or sets the discrete wavelet decomposition.
        /// </summary>
        public WaveletDecomposition WaveletDecomposition
        {
            get
            {
                return waveletDecomposition;
            }
            set
            {
                this.waveletDecomposition = value;
            }
        }
        #endregion

        #region Wavelet transform
        /// <summary>
        /// Forward wavelet transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            var packs = waveletDecomposition.Forward(A);
            int N = A.Length;
            var y = new float[N];
            int pos = 0;

            for (int i = 0; i < packs.Length; i++)
            {
                var b = packs[i];
                Buffer.BlockCopy(b, 0, y, pos * sizeof(float), b.Length * sizeof(float));
                pos += b.Length;
            }
            return y;
        }
        /// <summary>
        /// Backward wavelet transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length, L = Math.Min((int)Maths.Log2(N), waveletDecomposition.Levels);
            int aLen = N >> L; int pos = 0;
            var coeffs = new float[L + 1][];
            coeffs[0] = new float[aLen];
            Buffer.BlockCopy(B, 0, coeffs[0], 0, aLen * sizeof(float)); pos += aLen;

            for (int lev = L; lev >= 1; lev--)
            {
                int dLen = N >> lev; var d = new float[dLen];
                Buffer.BlockCopy(B, pos * sizeof(float), d, 0, dLen * sizeof(float)); pos += dLen; coeffs[L - lev + 1] = d;
            }

            return waveletDecomposition.Backward(coeffs);
        }
        /// <summary>
        /// Forward wavelet transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            var packs = waveletDecomposition.Forward(A); int L = (packs.Length - 1) / 3;
            int rows = packs[0].GetLength(0) << L;
            int cols = packs[0].GetLength(1) << L; var Y = new float[rows, cols];

            // place coarsest LL
            for (int r = 0; r < packs[0].GetLength(0); r++)
            {
                for (int c = 0; c < packs[0].GetLength(1); c++)
                {
                    Y[r, c] = packs[0][r, c];
                }
            }

            for (int lev = L; lev >= 1; lev--)
            {
                int h = rows >> (lev - 1), w = cols >> (lev - 1), hr = h >> 1, wc = w >> 1;
                var LH = packs[1 + 3 * (L - lev) + 0];
                var HL = packs[1 + 3 * (L - lev) + 1];
                var HH = packs[1 + 3 * (L - lev) + 2];

                for (int r = 0; r < hr; r++)
                {
                    for (int c = 0; c < wc; c++)
                    {
                        Y[r, wc + c] = LH[r, c];
                        Y[hr + r, c] = HL[r, c];
                        Y[hr + r, wc + c] = HH[r, c];
                    }
                }
            }

            return Y;
        }
        /// <summary>
        /// Backward wavelet transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            int R = B.GetLength(0), C = B.GetLength(1);
            int L = Math.Min(Math.Min((int)Maths.Log2(R), (int)Maths.Log2(C)), waveletDecomposition.Levels);
            var packs = new float[1 + 3 * L][,];
            int hL = R >> L, wL = C >> L; var LL = new float[hL, wL];

            for (int r = 0; r < hL; r++)
            {
                for (int c = 0; c < wL; c++)
                {
                    LL[r, c] = B[r, c]; packs[0] = LL;
                }
            }

            for (int lev = L; lev >= 1; lev--)
            {
                int h = R >> (lev - 1), w = C >> (lev - 1), hr = h >> 1, wc = w >> 1;
                var LH = new float[hr, wc];
                var HL = new float[hr, wc];
                var HH = new float[hr, wc];

                for (int r = 0; r < hr; r++)
                {
                    for (int c = 0; c < wc; c++)
                    {
                        LH[r, c] = B[r, wc + c];
                        HL[r, c] = B[hr + r, c];
                        HH[r, c] = B[hr + r, wc + c];
                    }
                }

                packs[1 + 3 * (L - lev) + 0] = LH;
                packs[1 + 3 * (L - lev) + 1] = HL;
                packs[1 + 3 * (L - lev) + 2] = HH;
            }

            return waveletDecomposition.Backward(packs);
        }
        /// <summary>
        /// Forward wavelet transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            var packs = waveletDecomposition.Forward(A);
            int N = A.Length;
            var y = new Complex32[N];
            int pos = 0;

            for (int i = 0; i < packs.Length; i++)
            {
                var b = packs[i];
                Array.Copy(b, 0, y, pos, b.Length);
                pos += b.Length;
            }

            return y;
        }
        /// <summary>
        /// Backward wavelet transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length, L = Math.Min((int)Maths.Log2(N), waveletDecomposition.Levels);
            int aLen = N >> L;
            int pos = 0;
            var coeffs = new Complex32[L + 1][];
            coeffs[0] = new Complex32[aLen];
            Array.Copy(B, pos, coeffs[0], 0, aLen);
            pos += aLen;

            for (int lev = L; lev >= 1; lev--)
            {
                int dLen = N >> lev;
                var d = new Complex32[dLen];
                Array.Copy(B, pos, d, 0, dLen);
                pos += dLen;
                coeffs[L - lev + 1] = d;
            }

            return waveletDecomposition.Backward(coeffs);
        }
        /// <summary>
        /// Forward wavelet transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            var packs = waveletDecomposition.Forward(A);
            int L = (packs.Length - 1) / 3;
            int rows = packs[0].GetLength(0) << L;
            int cols = packs[0].GetLength(1) << L;
            var Y = new Complex32[rows, cols];

            for (int r = 0; r < packs[0].GetLength(0); r++)
            {
                for (int c = 0; c < packs[0].GetLength(1); c++)
                {
                    Y[r, c] = packs[0][r, c];
                }
            }

            for (int lev = L; lev >= 1; lev--)
            {
                int h = rows >> (lev - 1), w = cols >> (lev - 1), hr = h >> 1, wc = w >> 1;
                var LH = packs[1 + 3 * (L - lev) + 0];
                var HL = packs[1 + 3 * (L - lev) + 1];
                var HH = packs[1 + 3 * (L - lev) + 2];

                for (int r = 0; r < hr; r++)
                {
                    for (int c = 0; c < wc; c++)
                    {
                        Y[r, wc + c] = LH[r, c];
                        Y[hr + r, c] = HL[r, c];
                        Y[hr + r, wc + c] = HH[r, c];
                    }
                }
            }

            return Y;
        }
        /// <summary>
        /// Backward wavelet transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int R = B.GetLength(0), C = B.GetLength(1);
            int L = Math.Min(Math.Min((int)Maths.Log2(R), (int)Maths.Log2(C)), waveletDecomposition.Levels);
            var packs = new Complex32[1 + 3 * L][,];
            int hL = R >> L, wL = C >> L;
            var LL = new Complex32[hL, wL];

            for (int r = 0; r < hL; r++)
            {
                for (int c = 0; c < wL; c++)
                {
                    LL[r, c] = B[r, c];
                }
            }

            packs[0] = LL;

            for (int lev = L; lev >= 1; lev--)
            {
                int h = R >> (lev - 1), w = C >> (lev - 1), hr = h >> 1, wc = w >> 1;
                var LH = new Complex32[hr, wc];
                var HL = new Complex32[hr, wc];
                var HH = new Complex32[hr, wc];

                for (int r = 0; r < hr; r++)
                {
                    for (int c = 0; c < wc; c++)
                    {
                        LH[r, c] = B[r, wc + c];
                        HL[r, c] = B[hr + r, c];
                        HH[r, c] = B[hr + r, wc + c];
                    }
                }

                packs[1 + 3 * (L - lev) + 0] = LH;
                packs[1 + 3 * (L - lev) + 1] = HL;
                packs[1 + 3 * (L - lev) + 2] = HH;
            }

            return waveletDecomposition.Backward(packs);
        }
        #endregion
    }
}