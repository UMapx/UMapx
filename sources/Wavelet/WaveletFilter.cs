using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the wavelet filter.
    /// </summary>
    [Serializable]
    public class WaveletFilter : IFilter
    {
        #region Private data
        WaveletTransform dwt;
        double factor = -1.0;
        double accuracy = 0.1;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the wavelet filter.
        /// </summary>
        /// <param name="dwt">Discrete wavelet transform</param>
        /// <param name="factor">Factor [-1, 1]</param>
        /// <param name="accuracy">Accuracy [0, 1]</param>
        public WaveletFilter(WaveletTransform dwt, double factor = -1.0, double accuracy = 0.1)
        {
            this.dwt = dwt;
            this.factor = factor;
            this.Accuracy = accuracy;
        }
        /// <summary>
        /// Gets or sets the discrete wavelet transform.
        /// </summary>
        public WaveletTransform DWT
        {
            get
            {
                return this.dwt;
            }
            set
            {
                this.dwt = value;
            }
        }
        /// <summary>
        /// Gets or sets the accuracy of the filter [0, 1].
        /// </summary>
        public double Accuracy
        {
            get
            {
                return this.accuracy;
            }
            set
            {
                this.accuracy = Maths.Double(value);
            }
        }
        /// <summary>
        /// Gets or sets the factor value [-1, 1].
        /// </summary>
        public double Factor
        {
            get
            {
                return this.factor;
            }
            set
            {
                this.factor = value;
            }
        }
        #endregion

        #region Public apply voids
        /// <summary>
        /// Implements a wavelet filter.
        /// </summary>
        /// <param name="data">Matrix</param>

        public void Apply(double[,] data)
        {
            // extend input
            int r0 = data.GetLength(0), c0 = data.GetLength(1);
            int delta = (int)(Math.Min(r0, c0) * accuracy);
            int r = WaveletFilter.GetLength(r0 + delta * 2, dwt.Levels);
            int c = WaveletFilter.GetLength(c0 + delta * 2, dwt.Levels);
            double[,] extd = Matrice.Extend(data, r, c);

            // wavelets
            double[,] wave;
            double alfa = 1 + factor;
            int maxl = WaveletFilter.GetMaxLevels(Math.Min(r, c), dwt.Levels);
            int powR = r >> maxl;
            int powC = c >> maxl;
            int j, k;

            // forward wavelet transform
            wave = dwt.Forward(extd);

            // do job
            for (j = 0; j < r; j++)
            {
                for (k = 0; k < c; k++)
                {
                    if (j < powR && k < powC)
                        // low-pass
                        extd[j, k] = wave[j, k];
                    else
                        // high-pass filter
                        extd[j, k] = wave[j, k] * alfa;
                }
            }

            // backward wavelet transform
            extd = dwt.Backward(extd);

            // cutend result
            int y0 = (r - r0) / 2;
            int x0 = (c - c0) / 2;
            extd = Matrice.Cut(extd, y0, x0, r0, c0);

            for (j = 0; j < r0; j++)
                for (k = 0; k < c0; k++)
                    data[j, k] = extd[j, k];

            return;
        }
        /// <summary>
        /// Implements a wavelet filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(double[] data)
        {
            // params
            int r0 = data.GetLength(0);
            int delta = (int)(r0 * accuracy);
            int r = WaveletFilter.GetLength(r0 + delta * 2, dwt.Levels);
            double[] extd = Matrice.Extend(data, r);

            // wavelets
            double[] wave;
            double alfa = 1 + factor;
            int maxl = WaveletFilter.GetMaxLevels(r, dwt.Levels);
            int powR = r >> maxl;
            int j;

            // forward wavelet transform
            wave = dwt.Forward(extd);

            // do job
            for (j = 0; j < r; j++)
            {
                if (j < powR)
                    // low-pass
                    extd[j] = wave[j];
                else
                    // high-pass filter
                    extd[j] = wave[j] * alfa;
            }

            // backward wavelet transform
            extd = dwt.Backward(extd);

            // cutend result
            int y0 = (r - r0) / 2;
            extd = Matrice.Cut(extd, y0, r0);

            for (j = 0; j < r0; j++)
                data[j] = extd[j];

            return;
        }
        /// <summary>
        /// Implements a wavelet filter.
        /// </summary>
        /// <param name="data">Matrix</param>

        public void Apply(Complex[,] data)
        {
            // extend input
            int r0 = data.GetLength(0), c0 = data.GetLength(1);
            int delta = (int)(Math.Min(r0, c0) * accuracy);
            int r = WaveletFilter.GetLength(r0 + delta * 2, dwt.Levels);
            int c = WaveletFilter.GetLength(c0 + delta * 2, dwt.Levels);
            Complex[,] extd = Matrice.Extend(data, r, c);

            // wavelets
            Complex[,] wave;
            double alfa = 1 + factor;
            int maxl = WaveletFilter.GetMaxLevels(Math.Min(r, c), dwt.Levels);
            int powR = r >> maxl;
            int powC = c >> maxl;
            int j, k;

            // forward wavelet transform
            wave = dwt.Forward(extd);

            // do job
            for (j = 0; j < r; j++)
            {
                for (k = 0; k < c; k++)
                {
                    if (j < powR && k < powC)
                        // low-pass
                        extd[j, k] = wave[j, k];
                    else
                        // high-pass filter
                        extd[j, k] = wave[j, k] * alfa;
                }
            }

            // backward wavelet transform
            extd = dwt.Backward(extd);

            // cutend result
            int y0 = (r - r0) / 2;
            int x0 = (c - c0) / 2;
            extd = Matrice.Cut(extd, y0, x0, r0, c0);

            for (j = 0; j < r0; j++)
                for (k = 0; k < c0; k++)
                    data[j, k] = extd[j, k];

            return;
        }
        /// <summary>
        /// Implements a wavelet filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(Complex[] data)
        {
            // params
            int r0 = data.GetLength(0);
            int delta = (int)(r0 * accuracy);
            int r = WaveletFilter.GetLength(r0 + delta * 2, dwt.Levels);
            Complex[] extd = Matrice.Extend(data, r);

            // wavelets
            Complex[] wave;
            double alfa = 1 + factor;
            int maxl = WaveletFilter.GetMaxLevels(r, dwt.Levels);
            int powR = r >> maxl;
            int j;

            // forward wavelet transform
            wave = dwt.Forward(extd);

            // do job
            for (j = 0; j < r; j++)
            {
                if (j < powR)
                    // low-pass
                    extd[j] = wave[j];
                else
                    // high-pass filter
                    extd[j] = wave[j] * alfa;
            }

            // backward wavelet transform
            extd = dwt.Backward(extd);

            // cutend result
            int y0 = (r - r0) / 2;
            extd = Matrice.Cut(extd, y0, r0);

            for (j = 0; j < r0; j++)
                data[j] = extd[j];

            return;
        }
        #endregion

        #region Static voids
        /// <summary>
        /// Returns the length value for transform.
        /// </summary>
        /// <param name="n">Length</param>
        /// <param name="levels">Number of levels</param>
        /// <returns>Length</returns>
        public static int GetLength(int n, int levels)
        {
            // params
            int log2 = GetMaxLevels(n, levels);
            int s = n, m, i;

            // do job
            for (i = 0; i < log2; i++)
            {
                m = Maths.Mod(s, 2);

                if (s >= 2)
                {
                    if (m != 0)
                        s = s + 1;
                    s = s / 2;
                }
            }

            return s * (int)Math.Pow(2, i);
        }
        /// <summary>
        /// Returns max levels of 2^K transform.
        /// </summary>
        /// <param name="n">Length</param>
        /// <param name="levels">Levels</param>
        /// <returns>New length</returns>
        private static int GetMaxLevels(int n, int levels)
        {
            return (int)Math.Min(Math.Log(n, 2), levels);
        }
        #endregion
    }
}
