using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the wavelet filter.
    /// <remarks>
    /// For the correct wavelet transform of a signal, it is necessary that its dimension be a power of 2.
    /// It is recommended to use Coiflets (C1, C2, C3, C4 and C5).
    /// </remarks>
    /// </summary>
    [Serializable]
    public class WaveletFilter : IFilter
    {
        #region Filter components
        /// <summary>
        /// Initializes the wavelet filter.
        /// </summary>
        /// <param name="waveletDecomposition">Discrete wavelet decomposition</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public WaveletFilter(WaveletDecomposition waveletDecomposition, float factor = -1.0f)
        {
            WaveletDecomposition = waveletDecomposition;
            Factor = factor;
        }
        /// <summary>
        /// Gets or sets the discrete wavelet transform.
        /// </summary>
        public WaveletDecomposition WaveletDecomposition { get; set; }
        /// <summary>
        /// Gets or sets the factor value [-1, 1].
        /// </summary>
        public float Factor { get; set; }
        #endregion

        #region Public apply voids
        /// <summary>
        /// Implements a wavelet filter.
        /// </summary>
        /// <param name="data">Matrix</param>

        public void Apply(float[,] data)
        {
            var B = WaveletDecomposition.Forward(data);
            int n, m;

            for (int i = 1; i < B.Length; i++)
            {
                var b = B[i];
                n = b.GetLength(0);
                m = b.GetLength(1);

                if (!Maths.IsEven(n)) n--; // ?
                if (!Maths.IsEven(m)) m--; // ?

                // visu_shrink
                var bb = Matrice.Reshape(b, b.Length).Abs();
                Array.Sort(bb);
                var median = Math.Sqrt(bb[bb.Length / 2]) * Math.Sqrt(Maths.Log2(n * m));

                for (int y = 0; y < n; y++)
                {
                    for (int x = 0; x < m; x++)
                    {
                        b[y, x] = Math.Abs(b[y, x]) > median ? b[y, x] : 0;
                    }
                }

                B[i] = b;
            }

            var output = WaveletDecomposition.Backward(B);
            float factor = 1 + Factor;
            n = data.GetLength(0);
            m = data.GetLength(1);

            for (int y = 0; y < n; y++)
            {
                for (int x = 0; x < m; x++)
                {
                    // pass filtering
                    data[y, x] = output[y, x] + factor * (data[y, x] - output[y, x]);
                }
            }
        }
        /// <summary>
        /// Implements a wavelet filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(float[] data)
        {
            var B = WaveletDecomposition.Forward(data);
            int n;

            for (int i = 1; i < B.Length; i++)
            {
                var b = B[i];
                n = b.GetLength(0);

                if (!Maths.IsEven(n)) n--; // ?

                // visu_shrink
                var bb = b.Abs();
                Array.Sort(bb);
                var median = Math.Sqrt(bb[bb.Length / 2]) * Math.Sqrt(Maths.Log2(n));

                for (int y = 0; y < n; y++)
                {
                    b[y] = Math.Abs(b[y]) > median ? b[y] : 0;
                }

                B[i] = b;
            }

            var output = WaveletDecomposition.Backward(B);
            float factor = 1 + Factor;
            n = data.GetLength(0);

            for (int y = 0; y < n; y++)
            {
                // pass filtering
                data[y] = output[y] + factor * (data[y] - output[y]);
            }
        }
        /// <summary>
        /// Implements a wavelet filter.
        /// </summary>
        /// <param name="data">Matrix</param>

        public void Apply(Complex32[,] data)
        {
            var B = WaveletDecomposition.Forward(data);
            int n, m;

            for (int i = 1; i < B.Length; i++)
            {
                var b = B[i];
                n = b.GetLength(0);
                m = b.GetLength(1);

                if (!Maths.IsEven(n)) n--; // ?
                if (!Maths.IsEven(m)) m--; // ?

                // visu_shrink
                var bb = Matrice.Reshape(b, b.Length).Abs();
                Array.Sort(bb);
                var median = Math.Sqrt(bb[bb.Length / 2]) * Math.Sqrt(Maths.Log2(n * m));

                for (int y = 0; y < n; y++)
                {
                    for (int x = 0; x < m; x++)
                    {
                        b[y, x] = Maths.Abs(b[y, x]) > median ? b[y, x] : 0;
                    }
                }

                B[i] = b;
            }

            var output = WaveletDecomposition.Backward(B);
            float factor = 1 + Factor;
            n = data.GetLength(0);
            m = data.GetLength(1);

            for (int y = 0; y < n; y++)
            {
                for (int x = 0; x < m; x++)
                {
                    // pass filtering
                    data[y, x] = output[y, x] + factor * (data[y, x] - output[y, x]);
                }
            }
        }
        /// <summary>
        /// Implements a wavelet filter.
        /// </summary>
        /// <param name="data">Array</param>

        public void Apply(Complex32[] data)
        {
            var B = WaveletDecomposition.Forward(data);
            int n;

            for (int i = 1; i < B.Length; i++)
            {
                var b = B[i];
                n = b.GetLength(0);

                if (!Maths.IsEven(n)) n--; // ?

                // visu_shrink
                var bb = b.Abs();
                Array.Sort(bb);
                var median = Math.Sqrt(bb[bb.Length / 2]) * Math.Sqrt(Maths.Log2(n));

                for (int y = 0; y < n; y++)
                {
                    b[y] = Maths.Abs(b[y]) > median ? b[y] : 0;
                }

                B[i] = b;
            }

            var output = WaveletDecomposition.Backward(B);
            float factor = 1 + Factor;
            n = data.GetLength(0);

            for (int y = 0; y < n; y++)
            {
                // pass filtering
                data[y] = output[y] + factor * (data[y] - output[y]);
            }
        }
        #endregion
    }
}
