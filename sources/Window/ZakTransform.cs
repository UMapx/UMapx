using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Zak transform.
    /// </summary>
    [Serializable]
    public class ZakTransform : IZakTransform
    {
        #region Private data
        private readonly FourierTransform DFT = new FourierTransform(false, Direction.Vertical);
        private int m;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Zak transform.
        /// </summary>
        /// <param name="m">Number of frequency shifts [4, N/2]</param>
        public ZakTransform(int m)
        {
            M = m;
        }
        /// <summary>
        /// Gets or sets number of frequency shifts [4, N/2].
        /// <remarks>
        /// Even number.
        /// </remarks>
        /// </summary>
        public int M
        {
            get
            {
                return m;
            }
            set
            {
                if (value <= 0) 
                    throw new ArgumentException("M must be positive");

                m = value;
            }
        }
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

            Complex32[] buffer = new Complex32[L];

            for (int m = 0; m < N; m++)
            {
                for (int k = 0; k < L; k++)
                {
                    int index = (m - M * k) % N;
                    if (index < 0) index += N;

                    buffer[k] = input[index];
                }

                buffer = DFT.Forward(buffer);

                for (int n = 0; n < L; n++)
                {
                    result[m, n] = buffer[n];
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

            Complex32[] buffer = new Complex32[L];

            for (int m = 0; m < N; m++)
            {
                for (int k = 0; k < L; k++)
                {
                    int index = (m - M * k) % N;
                    if (index < 0) index += N;

                    buffer[k] = input[index];
                }

                buffer = DFT.Forward(buffer);

                for (int n = 0; n < L; n++)
                {
                    result[m, n] = buffer[n];
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
            if (N % M != 0 || L != N / M) throw new ArgumentException("The dimensions of matrix do not correspond to the M parameter");

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
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public float[] Orthogonalize(float[] input)
        {
            int N = input.Length;
            if (N % M != 0) throw new ArgumentException("The length of the input must be a multiple of M");

            int L = N / M;
            int R = Maths.Gcd(M, L);           // number of Zak cosets to stitch
            int RL = R * L;                    // vertical length
            int delta = M / R;                 // time shift step between cosets
            const float eps = 1e-8f;

            // 1) Build G (RL x N) by circular shifts with step Δ = M/R
            var G = new Complex32[RL, N];
            for (int i = 0; i < RL; i++)
            {
                int shift = delta * i % N;
                for (int j = 0; j < N; j++)
                {
                    G[i, j] = input[Maths.Mod(j + shift, N)];
                }
            }

            // 2) Vertical Zak-DFT (length = RL)
            var Z = DFT.Forward(G); // Z has shape (RL x N)

            // 3) Block-wise normalization over R rows: i, i+L, i+2L, ...
            float w = R / (float)Math.Sqrt(M); // -> sum_r |Z|^2 becomes w^2 (matches R=2 case)
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    // accumulate energy over the R cosets at this (i,j)
                    float energy = 0f;
                    for (int r = 0; r < R; r++)
                    {
                        var z = Z[i + r * L, j];
                        energy += z.Real * z.Real + z.Imag * z.Imag;
                    }
                    float phi = w / (float)Math.Sqrt(energy + eps);

                    // scale all R components by the same φ
                    for (int r = 0; r < R; r++)
                        Z[i + r * L, j] *= phi;
                }
            }

            // 4) Take DC over the vertical axis (average across RL rows)
            var outR = new float[N];
            for (int j = 0; j < N; j++)
            {
                Complex32 sum = Complex32.Zero;
                for (int i = 0; i < RL; i++)
                    sum += Z[i, j];
                outR[j] = (sum / RL).Real;
            }
            return outR;
        }
        /// <summary>
        /// Zak orthogonalization.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Orthogonalize(Complex32[] input)
        {
            int N = input.Length;
            if (N % M != 0) throw new ArgumentException("The length of the input must be a multiple of M");

            int L = N / M;
            int R = Maths.Gcd(M, L);           // number of Zak cosets to stitch
            int RL = R * L;                    // vertical length
            int delta = M / R;                 // time shift step between cosets
            const float eps = 1e-8f;

            // 1) Build G (RL x N) by circular shifts with step Δ = M/R
            var G = new Complex32[RL, N];
            for (int i = 0; i < RL; i++)
            {
                int shift = delta * i % N;
                for (int j = 0; j < N; j++)
                {
                    G[i, j] = input[Maths.Mod(j + shift, N)];
                }
            }

            // 2) Vertical Zak-DFT (length = RL)
            var Z = DFT.Forward(G); // Z has shape (RL x N)

            // 3) Block-wise normalization over R rows: i, i+L, i+2L, ...
            float w = R / (float)Math.Sqrt(M); // -> sum_r |Z|^2 becomes w^2 (matches R=2 case)
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    // accumulate energy over the R cosets at this (i,j)
                    float energy = 0f;
                    for (int r = 0; r < R; r++)
                    {
                        var z = Z[i + r * L, j];
                        energy += z.Real * z.Real + z.Imag * z.Imag;
                    }
                    float phi = w / (float)Math.Sqrt(energy + eps);
                    for (int r = 0; r < R; r++)
                        Z[i + r * L, j] *= phi;
                }
            }

            // 4) Take DC over the vertical axis (average across RL rows)
            var outC = new Complex32[N];
            for (int j = 0; j < N; j++)
            {
                Complex32 sum = Complex32.Zero;
                for (int i = 0; i < RL; i++)
                    sum += Z[i, j];
                outC[j] = sum / RL;
            }
            return outC;
        }
        #endregion
    }
}
