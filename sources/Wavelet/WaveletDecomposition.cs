using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines a discrete wavelet decomposition.
    /// <remarks>
    /// For the correct wavelet transform of a signal, it is necessary that its dimension be a power of 2.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_wavelet_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class WaveletDecomposition : IPyramidTransform
    {
        #region Initialize
        /// <summary>
        /// Initializes a discrete wavelet decomposition.
        /// </summary>
        /// <param name="waveletTransform">Discrete wavelet transform</param>
        public WaveletDecomposition(WaveletTransform waveletTransform)
        {
            WaveletTransform = waveletTransform;
        }
        /// <summary>
        /// Gets or sets the discrete wavelet.
        /// </summary>
        public WaveletTransform WaveletTransform { get; set; }
        #endregion

        #region Wavelet transform
        /// <summary>
        /// Forward wavelet decomposition.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[][] Forward(float[] A)
        {
            // params
            int length = A.Length;
            int nLevels = (int)Math.Min(Maths.Log2(length), WaveletTransform.Levels);
            float[] B = WaveletTransform.Forward(A);
            float[][] C = new float[nLevels + 1][];

            // forward multi-scale wavelet decomposition
            for (int i = 0, k = 0; i <= nLevels; i++)
            {
                var bound = length >> (nLevels - i);
                var count = bound - k;

                C[i] = new float[count];

                for (int j = 0; j < count; j++)
                {
                    C[i][j] = B[j + k];
                }

                k = bound;
            }

            return C;
        }
        /// <summary>
        /// Backward wavelet decomposition.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[][] B)
        {
            // params
            int nLevels = B.Length;
            float[] A = new float[] { };

            // backward multi-scale wavelet decomposition
            for (int i = 0; i < nLevels; i++)
            {
                A = Matrice.Merge(A, B[i]);
            }

            return WaveletTransform.Backward(A);
        }
        /// <summary>
        /// Forward wavelet decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[][,] Forward(float[,] A)
        {
            // params
            int N = A.GetLength(0);
            int M = A.GetLength(1);
            int nLevels = (int)Math.Min(Math.Min(Maths.Log2(N), WaveletTransform.Levels), M);
            float[,] B = WaveletTransform.Forward(A);
            float[][,] C = new float[nLevels + 1][,];

            // forward multi-scale wavelet decomposition
            for (int i = 0, k1 = 0, k2 = 0; i <= nLevels; i++)
            {
                var bound1 = N >> (nLevels - i);
                var bound2 = M >> (nLevels - i);
                var count1 = bound1 - k1;
                var count2 = bound2 - k2;

                C[i] = new float[count1, count2];

                for (int y = 0; y < count1; y++)
                {
                    for (int x = 0; x < count2; x++)
                    {
                        C[i][y, x] = B[y + k1, x + k2];
                        B[y + k1, x + k2] = 0;
                    }
                }

                //k1 = bound1; // not actual for 2D decomposition
                //k2 = bound2;
            }

            return C;
        }
        /// <summary>
        /// Backward wavelet decomposition.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[][,] B)
        {
            // params
            int nLevels = B.Length;
            float[,] A = B[nLevels - 1];

            // backward multi-scale wavelet decomposition
            for (int i = nLevels - 2; i >= 0; i--)
            {
                var b = B[i];
                var col = b.GetLength(0);
                var row = b.GetLength(1);

                for (int y = 0; y < col; y++)
                {
                    for (int x = 0; x < row; x++)
                    {
                        A[y, x] = b[y, x];
                    }
                }
            }

            return WaveletTransform.Backward(A);
        }
        /// <summary>
        /// Forward wavelet decomposition.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[][] Forward(Complex32[] A)
        {
            // params
            int length = A.Length;
            int nLevels = (int)Math.Min(Maths.Log2(length), WaveletTransform.Levels);
            Complex32[] B = WaveletTransform.Forward(A);
            Complex32[][] C = new Complex32[nLevels + 1][];

            // forward multi-scale wavelet decomposition
            for (int i = 0, k = 0; i <= nLevels; i++)
            {
                var bound = length >> (nLevels - i);
                var count = bound - k;

                C[i] = new Complex32[count];

                for (int j = 0; j < count; j++)
                {
                    C[i][j] = B[j + k];
                }

                k = bound;
            }

            return C;
        }
        /// <summary>
        /// Backward wavelet decomposition.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[][] B)
        {
            // params
            int nLevels = B.Length;
            Complex32[] A = new Complex32[] { };

            // backward multi-scale wavelet decomposition
            for (int i = 0; i < nLevels; i++)
            {
                A = Matrice.Merge(A, B[i]);
            }

            return WaveletTransform.Backward(A);
        }
        /// <summary>
        /// Forward wavelet decomposition.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[][,] Forward(Complex32[,] A)
        {
            // params
            int N = A.GetLength(0);
            int M = A.GetLength(1);
            int nLevels = (int)Math.Min(Math.Min(Maths.Log2(N), WaveletTransform.Levels), M);
            Complex32[,] B = WaveletTransform.Forward(A);
            Complex32[][,] C = new Complex32[nLevels + 1][,];

            // forward multi-scale wavelet decomposition
            for (int i = 0, k1 = 0, k2 = 0; i <= nLevels; i++)
            {
                var bound1 = N >> (nLevels - i);
                var bound2 = M >> (nLevels - i);
                var count1 = bound1 - k1;
                var count2 = bound2 - k2;

                C[i] = new Complex32[count1, count2];

                for (int y = 0; y < count1; y++)
                {
                    for (int x = 0; x < count2; x++)
                    {
                        C[i][y, x] = B[y + k1, x + k2];
                        B[y + k1, x + k2] = 0;
                    }
                }

                //k1 = bound1; // not actual for 2D decomposition
                //k2 = bound2;
            }

            return C;
        }
        /// <summary>
        /// Backward wavelet decomposition.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[][,] B)
        {
            // params
            int nLevels = B.Length;
            Complex32[,] A = B[nLevels - 1];

            // backward multi-scale wavelet decomposition
            for (int i = nLevels - 2; i >= 0; i--)
            {
                var b = B[i];
                var col = b.GetLength(0);
                var row = b.GetLength(1);

                for (int y = 0; y < col; y++)
                {
                    for (int x = 0; x < row; x++)
                    {
                        A[y, x] = b[y, x];
                    }
                }
            }

            return WaveletTransform.Backward(A);
        }
        #endregion
    }
}
