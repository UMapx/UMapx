using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines a discrete wavelet transform.
    /// <remarks>
    /// For the correct wavelet transform of a signal, it is necessary that its dimension be a power of 2.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_wavelet_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class WaveletTransform : IWaveletTransform, ITransform
    {
        #region Private data
        private float[] lp;        // Low-Pass filter,
        private float[] hp;        // High-Pass filer,
        private float[] ilp;       // Inverse Low-Pass filter,
        private float[] ihp;       // Inverse High-Pass filter,
        private bool normalized;    // Normalized transform or not,
        private int levels;         // Number of levels.
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes a discrete wavelet transform.
        /// </summary>
        /// <param name="wavelet">Discrete wavelet</param>
        /// <param name="levels">Number of levels</param>
        /// <param name="normalized">Normalized transform or not</param>
        public WaveletTransform(WaveletPack wavelet, int levels = 1, bool normalized = true)
        {
            Wavelet = wavelet; Levels = levels; Normalized = normalized;
        }
        /// <summary>
        /// Gets or sets the number of transform levels.
        /// </summary>
        public int Levels
        {
            get
            {
                return this.levels;
            }
            set
            {
                if (value < 1)
                    throw new Exception("Number of levels cannot be less than 1");

                this.levels = value;
            }
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
        /// <summary>
        /// Gets or sets the discrete wavelet.
        /// </summary>
        public WaveletPack Wavelet
        {
            get
            {
                return new WaveletPack(lp, hp, ilp, ihp);
            }
            set
            {
                this.lp = value.LowPass;
                this.hp = value.HighPass;
                this.ilp = value.ILowPass;
                this.ihp = value.IHighPass;
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
            // params
            int nLevels = (int)Math.Min(Maths.Log2(A.Length), this.levels);

            // forward multi-scale wavelet transform
            for (int i = 0; i < nLevels; i++)
            {
                A = this.dwt(A, i);
            }

            return A;
        }
        /// <summary>
        /// Backward wavelet transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            // params
            int nLevels = (int)Math.Min(Maths.Log2(B.Length), this.levels);

            // backward multi-scale wavelet transform
            for (int i = nLevels; i > 0; i--)
            {
                B = this.idwt(B, i);
            }

            return B;
        }
        /// <summary>
        /// Forward wavelet transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            // params
            int Bound1, Bound2, i, j;
            int DataLen1 = A.GetLength(0);
            int DataLen2 = A.GetLength(1);
            float[,] output = (float[,])A.Clone();
            float[] buff2 = new float[DataLen2];
            float[] buff1 = new float[DataLen1];
            int nLevels = (int)Math.Min(Math.Min(Maths.Log2(DataLen1), this.levels), DataLen2);

            // do job
            for (int lev = 0; lev < nLevels; lev++)
            {
                Bound1 = DataLen1 >> lev;
                Bound2 = DataLen2 >> lev;

                if (!Maths.IsEven(Bound1) && Bound1 < DataLen1)
                    Bound1--;
                if (!Maths.IsEven(Bound2) && Bound2 < DataLen2)
                    Bound2--;

                for (i = 0; i < Bound1; i++)
                {
                    for (j = 0; j < Bound2; j++) buff2[j] = output[i, j];
                    buff2 = this.dwt(buff2, lev);
                    for (j = 0; j < Bound2; j++) output[i, j] = buff2[j];
                }

                for (j = 0; j < Bound2; j++)
                {
                    for (i = 0; i < Bound1; i++) buff1[i] = output[i, j];
                    buff1 = this.dwt(buff1, lev);
                    for (i = 0; i < Bound1; i++) output[i, j] = buff1[i];
                }
            }

            return output;
        }
        /// <summary>
        /// Backward wavelet transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            // params
            int Bound1, Bound2, i, j;
            int DataLen1 = B.GetLength(0);
            int DataLen2 = B.GetLength(1);
            float[,] output = (float[,])B.Clone();
            float[] buff1 = new float[DataLen1];
            float[] buff2 = new float[DataLen2];
            int nLevels = (int)Math.Min(Math.Min(Maths.Log2(DataLen1), this.levels), DataLen2);

            // do job
            for (int lev = nLevels; lev > 0; lev--)
            {
                Bound1 = DataLen1 >> lev;
                Bound2 = DataLen2 >> lev;

                for (i = 0; i < Bound1 << 1; i++)
                {
                    for (j = 0; j < Bound2 << 1; j++) buff2[j] = output[i, j];
                    buff2 = this.idwt(buff2, lev);
                    for (j = 0; j < Bound2 << 1; j++) output[i, j] = buff2[j];
                }

                for (j = 0; j < Bound2 << 1; j++)
                {
                    for (i = 0; i < Bound1 << 1; i++) buff1[i] = output[i, j];
                    buff1 = this.idwt(buff1, lev);
                    for (i = 0; i < Bound1 << 1; i++) output[i, j] = buff1[i];
                }
            }

            return output;
        }
        /// <summary>
        /// Forward wavelet transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            // params
            int nLevels = (int)Math.Min(Maths.Log2(A.Length), this.levels);

            // forward multi-scale wavelet transform
            for (int i = 0; i < nLevels; i++)
            {
                A = this.dwt(A, i);
            }

            return A;
        }
        /// <summary>
        /// Backward wavelet transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            // params
            int nLevels = (int)Math.Min(Maths.Log2(B.Length), this.levels);

            // backward multi-scale wavelet transform
            for (int i = nLevels; i > 0; i--)
            {
                B = this.idwt(B, i);
            }

            return B;
        }
        /// <summary>
        /// Forward wavelet transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            // params
            int Bound1, Bound2, i, j;
            int DataLen1 = A.GetLength(0);
            int DataLen2 = A.GetLength(1);
            Complex32[,] output = (Complex32[,])A.Clone();
            Complex32[] buff2 = new Complex32[DataLen2];
            Complex32[] buff1 = new Complex32[DataLen1];
            int nLevels = (int)Math.Min(Math.Min(Maths.Log2(DataLen1), this.levels), DataLen2);

            // do job
            for (int lev = 0; lev < nLevels; lev++)
            {
                Bound1 = DataLen1 >> lev;
                Bound2 = DataLen2 >> lev;

                if (!Maths.IsEven(Bound1) && Bound1 < DataLen1)
                    Bound1--;
                if (!Maths.IsEven(Bound2) && Bound2 < DataLen2)
                    Bound2--;

                for (i = 0; i < Bound1; i++)
                {
                    for (j = 0; j < Bound2; j++) buff2[j] = output[i, j];
                    buff2 = this.dwt(buff2, lev);
                    for (j = 0; j < Bound2; j++) output[i, j] = buff2[j];
                }

                for (j = 0; j < Bound2; j++)
                {
                    for (i = 0; i < Bound1; i++) buff1[i] = output[i, j];
                    buff1 = this.dwt(buff1, lev);
                    for (i = 0; i < Bound1; i++) output[i, j] = buff1[i];
                }
            }

            return output;
        }
        /// <summary>
        /// Backward wavelet transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            // params
            int Bound1, Bound2, i, j;
            int DataLen1 = B.GetLength(0);
            int DataLen2 = B.GetLength(1);
            Complex32[,] output = (Complex32[,])B.Clone();
            Complex32[] buff1 = new Complex32[DataLen1];
            Complex32[] buff2 = new Complex32[DataLen2];
            int nLevels = (int)Math.Min(Math.Min(Maths.Log2(DataLen1), this.levels), DataLen2);

            // do job
            for (int lev = nLevels; lev > 0; lev--)
            {
                Bound1 = DataLen1 >> lev;
                Bound2 = DataLen2 >> lev;

                for (i = 0; i < Bound1 << 1; i++)
                {
                    for (j = 0; j < Bound2 << 1; j++) buff2[j] = output[i, j];
                    buff2 = this.idwt(buff2, lev);
                    for (j = 0; j < Bound2 << 1; j++) output[i, j] = buff2[j];
                }

                for (j = 0; j < Bound2 << 1; j++)
                {
                    for (i = 0; i < Bound1 << 1; i++) buff1[i] = output[i, j];
                    buff1 = this.idwt(buff1, lev);
                    for (i = 0; i < Bound1 << 1; i++) output[i, j] = buff1[i];
                }
            }

            return output;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Forward discrete wavelet transform.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="level">Current level of transform</param>
        /// <returns>Output data</returns>
        private float[] dwt(float[] input, int level)
        {
            // params
            int length = input.Length;
            float[] output = new float[length];
            int Bound = length >> level;

            // odd element
            if (!Maths.IsEven(Bound))
            {
                Bound--;
            }

            int lpLen = this.lp.Length;
            int hpLen = this.hp.Length;
            int lpStart = -((lpLen >> 1) - 1);
            int hpStart = -((hpLen >> 1) - 1);
            Array.Copy(input, Bound, output, Bound, length - Bound);
            float a = 0;
            float b = 0;
            int h = Bound >> 1;
            int c, i, j, r, k;

            // do job
            for (i = 0, r = 0; i < Bound; i += 2, r++)
            {
                // low-pass filter
                for (j = lpStart, k = 0; k < lpLen; j++, k++)
                {
                    if (j < 0 || j >= Bound)
                        c = (j % Bound + Bound) % Bound;
                    else
                        c = j;
                    a += this.lp[k] * input[c];
                }

                // high-pass filter
                for (j = hpStart, k = 0; k < hpLen; j++, k++)
                {
                    if (j < 0 || j >= Bound)
                        c = (j % Bound + Bound) % Bound;
                    else
                        c = j;
                    b += this.hp[k] * input[c];
                }
                lpStart += 2;
                hpStart += 2;

                if (normalized)
                {
                    output[r] = a / Maths.Sqrt2;
                    output[r + h] = b / Maths.Sqrt2;
                }
                else
                {
                    output[r] = a;
                    output[r + h] = b;
                }

                a = 0;
                b = 0;
            }

            return output;
        }
        /// <summary>
        /// Backward discrete wavelet transform.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="level">Current level of transform</param>
        /// <returns>Output data</returns>
        private float[] idwt(float[] input, int level)
        {
            // params
            int length = input.Length;
            float[] output = (float[])input.Clone();
            int Bound = length >> level;
            int h = Bound << 1;
            int lpLen = this.ilp.Length;
            int hpLen = this.ihp.Length;
            int lpStart = -((lpLen >> 1) - 1);
            int hpStart = -((hpLen >> 1) - 1);
            float[] Low = new float[h];
            float[] Hig = new float[h];
            float s = 0;
            int c, i, j, k;

            // redim
            for (i = 0, j = 0; i < h; i += 2, j++)
            {
                Low[i] = 0;
                Hig[i] = 0;
                Low[i + 1] = input[j];
                Hig[i + 1] = input[Bound + j];
            }

            // do job
            for (i = 0; i < h; i++)
            {
                // low-pass filter
                for (j = lpStart, k = 0; k < lpLen; j++, k++)
                {
                    if (j < 0 || j >= h)
                        c = (j % h + h) % h;
                    else
                        c = j;
                    s += this.ilp[k] * Low[c];
                }

                // high-pass filter
                for (j = hpStart, k = 0; k < hpLen; j++, k++)
                {
                    if (j < 0 || j >= h)
                        c = (j % h + h) % h;
                    else
                        c = j;
                    s += this.ihp[k] * Hig[c];
                }

                lpStart += 1;
                hpStart += 1;
                output[i] = (normalized) ? s * Maths.Sqrt2 : s;
                s = 0;
            }

            return output;
        }

        /// <summary>
        /// Forward discrete wavelet transform.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="level">Current level of transform</param>
        /// <returns>Output data</returns>
        private Complex32[] dwt(Complex32[] input, int level)
        {
            // params
            int length = input.Length;
            Complex32[] output = new Complex32[length];
            int Bound = length >> level;

            // odd element
            if (!Maths.IsEven(Bound))
            {
                Bound--;
            }

            int lpLen = this.lp.Length;
            int hpLen = this.hp.Length;
            int lpStart = -((lpLen >> 1) - 1);
            int hpStart = -((hpLen >> 1) - 1);
            Array.Copy(input, Bound, output, Bound, length - Bound);
            Complex32 a = 0;
            Complex32 b = 0;
            int h = Bound >> 1;
            int c, i, j, r, k;

            // do job
            for (i = 0, r = 0; i < Bound; i += 2, r++)
            {
                // low-pass filter
                for (j = lpStart, k = 0; k < lpLen; j++, k++)
                {
                    if (j < 0 || j >= Bound)
                        c = (j % Bound + Bound) % Bound;
                    else
                        c = j;
                    a += this.lp[k] * input[c];
                }

                // high-pass filter
                for (j = hpStart, k = 0; k < hpLen; j++, k++)
                {
                    if (j < 0 || j >= Bound)
                        c = (j % Bound + Bound) % Bound;
                    else
                        c = j;
                    b += this.hp[k] * input[c];
                }
                lpStart += 2;
                hpStart += 2;

                if (normalized)
                {
                    output[r] = a / Maths.Sqrt2;
                    output[r + h] = b / Maths.Sqrt2;
                }
                else
                {
                    output[r] = a;
                    output[r + h] = b;
                }

                a = 0;
                b = 0;
            }

            return output;
        }
        /// <summary>
        /// Backward discrete wavelet transform.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="level">Current level of transform</param>
        /// <returns>Output data</returns>
        private Complex32[] idwt(Complex32[] input, int level)
        {
            // params
            int length = input.Length;
            Complex32[] output = (Complex32[])input.Clone();
            int Bound = length >> level;
            int h = Bound << 1;
            int lpLen = this.ilp.Length;
            int hpLen = this.ihp.Length;
            int lpStart = -((lpLen >> 1) - 1);
            int hpStart = -((hpLen >> 1) - 1);
            Complex32[] Low = new Complex32[h];
            Complex32[] Hig = new Complex32[h];
            Complex32 s = 0;
            int c, i, j, k;

            // redim
            for (i = 0, j = 0; i < h; i += 2, j++)
            {
                Low[i] = 0;
                Hig[i] = 0;
                Low[i + 1] = input[j];
                Hig[i + 1] = input[Bound + j];
            }

            // do job
            for (i = 0; i < h; i++)
            {
                // low-pass filter
                for (j = lpStart, k = 0; k < lpLen; j++, k++)
                {
                    if (j < 0 || j >= h)
                        c = (j % h + h) % h;
                    else
                        c = j;
                    s += this.ilp[k] * Low[c];
                }

                // high-pass filter
                for (j = hpStart, k = 0; k < hpLen; j++, k++)
                {
                    if (j < 0 || j >= h)
                        c = (j % h + h) % h;
                    else
                        c = j;
                    s += this.ihp[k] * Hig[c];
                }

                lpStart += 1;
                hpStart += 1;
                output[i] = (normalized) ? s * Maths.Sqrt2 : s;
                s = 0;
            }

            return output;
        }
        #endregion
    }
}
