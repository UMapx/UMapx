using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Zak transform.
    /// </summary>
    [Serializable]
    public class ZakTransform
    {
        #region Private data
        private readonly FourierTransform DFT = new FourierTransform(false, Direction.Vertical);
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Zak transform.
        /// </summary>
        /// <param name="m">Number of frequency shifts [4, N]</param>
        public ZakTransform(int m)
        {
            M = m;
        }
        /// <summary>
        /// Gets or sets number of frequency shifts [4, N].
        /// <remarks>
        /// Even number.
        /// </remarks>
        /// </summary>
        public int M { get; set; }
        #endregion

        #region Zak transform
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Matrix</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        public Complex32[,] Forward(Complex32[] input)
        {
            int N = input.Length;
            if (N % M != 0) throw new ArgumentException("The length of the input must be a multiple of M");

            int L = N / M;
            var result = new Complex32[N, L];

            for (int m = 0; m < N; m++)
            {
                for (int n = 0; n < L; n++)
                {
                    Complex32 sum = Complex32.Zero;

                    for (int k = 0; k < L; k++)
                    {
                        int index = (m - M * k) % N;
                        if (index < 0) index += N;

                        float angle = 2 * Maths.Pi * n * k / L;
                        Complex32 w = Maths.Exp(new Complex32(0, angle)); // e^{i*angle}

                        sum += input[index] * w;
                    }

                    result[m, n] = sum;
                }
            }

            return result;
        }
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Matrix</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        public Complex32[,] Forward(float[] input)
        {
            int N = input.Length;
            if (N % M != 0) throw new ArgumentException("The length of the input must be a multiple of M");

            int L = N / M;
            var result = new Complex32[N, L];

            for (int m = 0; m < N; m++)
            {
                for (int n = 0; n < L; n++)
                {
                    Complex32 sum = Complex32.Zero;

                    for (int k = 0; k < L; k++)
                    {
                        int index = (m - M * k) % N;
                        if (index < 0) index += N;

                        float angle = 2 * Maths.Pi * n * k / L;
                        Complex32 w = Maths.Exp(new Complex32(0, angle)); // e^{i*angle}

                        sum += input[index] * w;
                    }

                    result[m, n] = sum;
                }
            }

            return result;
        }

        /// <summary>
        /// Backward Zak transform.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <returns>Array</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        public Complex32[] Backward(Complex32[,] matrix)
        {
            int N = matrix.GetLength(0);
            int L = matrix.GetLength(1);

            if (N % M != 0 || L != N / M)
                throw new ArgumentException("The dimensions of matrix do not correspond to the M parameter");

            var output = new Complex32[N];

            for (int m = 0; m < N; m++)
            {
                Complex32 sum = Complex32.Zero;

                for (int n = 0; n < L; n++)
                {
                    sum += matrix[m, n];
                }

                output[m] = sum / L;
            }

            return output;
        }
        #endregion

        #region Zak orthogonalization
        /// <summary>
        /// Zak orthogonalization.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Orthogonalize(float[] A)
        {
            // Fast shaping orthogonalization algorithm
            // WH functions using a discrete Zak transform.
            // V.P. Volchkov, D.A. Petrov and V.M. Asiryan.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = A.Length;
            float[] vort = new float[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex32[,] G = new Complex32[L2, N];
            Complex32[,] Z;

            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = A[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            Z = DFT.Forward(G);

            float w = 2 / (float)Math.Sqrt(M);
            float even, odd, phi;
            Complex32 z1, z2;

            for (i = 0; i < L; i++)
            {
                for (j = 0; j < N; j++)
                {
                    z1 = Z[i, j];
                    z2 = Z[L + i, j];

                    even = (float)Math.Pow(z1.Abs, 2);
                    odd = (float)Math.Pow(z2.Abs, 2);
                    phi = w / (float)Math.Sqrt(even + odd);

                    Z[i, j] = z1 * phi;
                    Z[L + i, j] = z2 * phi;
                }
            }

            Complex32 sum;
            for (i = 0; i < N; i++)
            {
                sum = 0;
                for (j = 0; j < L2; j++)
                {
                    sum += Z[j, i];
                }
                vort[i] = (sum / L2).Real;
            }

            return vort;
        }
        /// <summary>
        /// Zak orthogonalization.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Orthogonalize(Complex32[] A)
        {
            // Fast shaping orthogonalization algorithm
            // WH functions using a discrete Zak transform.
            // V.P. Volchkov, D.A. Petrov and V.M. Asiryan.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = A.Length;
            Complex32[] vort = new Complex32[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex32[,] G = new Complex32[L2, N];
            Complex32[,] Z;

            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = A[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            Z = DFT.Forward(G);

            float w = 2 / (float)Math.Sqrt(M);
            float even, odd, phi;
            Complex32 z1, z2;

            for (i = 0; i < L; i++)
            {
                for (j = 0; j < N; j++)
                {
                    z1 = Z[i, j];
                    z2 = Z[L + i, j];

                    even = (float)Math.Pow(z1.Abs, 2);
                    odd = (float)Math.Pow(z2.Abs, 2);
                    phi = w / (float)Math.Sqrt(even + odd);

                    Z[i, j] = z1 * phi;
                    Z[L + i, j] = z2 * phi;
                }
            }

            Complex32 sum;
            for (i = 0; i < N; i++)
            {
                sum = 0;
                for (j = 0; j < L2; j++)
                {
                    sum += Z[j, i];
                }
                vort[i] = sum / L2;
            }

            return vort;
        }
        #endregion
    }
}
