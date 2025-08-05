using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines fast Weyl-Heisenberg transform.
    /// <remarks>
    /// The class represents a computationally efficient implementation of one-dimensional and two-dimensional discrete orthogonal
    /// Weyl-Heisenberg transforms.
    /// More information can be found on the website:
    /// https://ieeexplore.ieee.org/document/9117707/
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastWeylHeisenbergTransform : IWindowTransform, ITransform
    {
        #region Private data
        private static readonly FastFourierTransform FFT = new FastFourierTransform(false, Direction.Horizontal);
        private IWindow window;
        private int m;
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes fast Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N]</param>
        /// <param name="direction">Processing direction</param>
        public FastWeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical)
        {
            Window = window; M = m; Direction = direction;
        }
        /// <summary>
        /// Gets or sets number of frequency shifts [4, N].
        /// <remarks>
        /// Even number.
        /// </remarks>
        /// </summary>
        public int M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value <= 2 || !Maths.IsEven(value))
                    throw new Exception("Invalid argument value");

                this.m = value;
            }
        }
        /// <summary>
        /// Gets or sets the processing direction.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Gets or sets the window function.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            float[] g0 = WeylHeisenbergTransform.GetPacket(this.window, A.Length);
            FastZakTransform zakTransform = new FastZakTransform(m);
            return FastWeylHeisenbergTransform.WHT(A, zakTransform.Orthogonalize(g0), m);
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            float[] g0 = WeylHeisenbergTransform.GetPacket(this.window, B.Length);
            FastZakTransform zakTransform = new FastZakTransform(m);
            return FastWeylHeisenbergTransform.IWHT(B, zakTransform.Orthogonalize(g0), m);
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            Complex32[,] B = (Complex32[,])A.Clone();
            int N = B.GetLength(0), M = B.GetLength(1);

            FastZakTransform zakTransform = new FastZakTransform(m);
            float[] g0 = zakTransform.Orthogonalize(WeylHeisenbergTransform.GetPacket(this.window, N));
            float[] g1 = zakTransform.Orthogonalize(WeylHeisenbergTransform.GetPacket(this.window, M));

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = WHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = WHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = WHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = WHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            Complex32[,] A = (Complex32[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            FastZakTransform zakTransform = new FastZakTransform(m);
            float[] g0 = zakTransform.Orthogonalize(WeylHeisenbergTransform.GetPacket(this.window, N));
            float[] g1 = zakTransform.Orthogonalize(WeylHeisenbergTransform.GetPacket(this.window, M));

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = IWHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = IWHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = IWHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = IWHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Public static components
        /// <summary>
        /// Fast forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
        public static Complex32[] WHT(Complex32[] input, float[] g0, int M)
        {
            // The function implements a fast Weil-Heisenberg direct transformation algorithm,
            // stated in the following articles:
            // A. Vahlin, "EFFICIENT ALGORITHMS FOR MODULATION AND DEMODULATION IN OFDM-SYSTEMS" [1].
            // V.M. Asiryan, V.P. Volchkov, "EFFECTIVE IMPLEMENTATION OF THE DIRECT TRANSFORMATION OF WEIL-HEISENBERG" [2].
            // The algorithm is computationally efficient for large M.

            int N = input.Length, L = N / M, M2 = M / 2, M4 = M2 / 2;
            Complex32[] output = new Complex32[N];
            Complex32[] exp = FastWeylHeisenbergTransform.GetRotation(M);

            Complex32[,] s0 = new Complex32[M, L];
            Complex32[,] a0 = new Complex32[M, L];
            Complex32[,] b0 = new Complex32[M, L];
            Complex32[,] A0 = new Complex32[L, M];
            Complex32[,] B0 = new Complex32[L, M];
            Complex32[,] A1 = new Complex32[L, M2];
            Complex32[,] B1 = new Complex32[L, M2];
            Complex32 c1re, c2re;
            Complex32 c1im, c2im;
            int k, i, j, u, n, m, l;

            for (m = 0; m < M; m++)
            {
                for (n = 0; n < L; n++)
                {
                    u = n * M;
                    i = Maths.Mod(m + M4 + u, N);
                    j = Maths.Mod(m - M4 + u, N);
                    k = Maths.Mod(-m - M4 - u, N);

                    s0[m, n] = input[k];
                    a0[m, n] = g0[i];
                    b0[m, n] = g0[j];
                }
            }

            for (l = 0; l < L; l++)
            {
                for (n = 0; n < L; n++)
                {
                    k = Maths.Mod(n - l, L);

                    for (m = 0; m < M; m++)
                    {
                        A0[l, m] += a0[m, n] * s0[m, k];
                        B0[l, m] += b0[m, n] * s0[m, k];
                    }
                }
            }

            Complex32 x, y, z, w;

            for (l = 0; l < L; l++)
            {
                for (m = 0; m < M2; m++)
                {
                    x = A0[l, m];
                    y = A0[l, m + M2];
                    z = A0[l, M2 - m].Conjugate;
                    w = A0[l, Maths.Mod(M - m, M)].Conjugate;

                    c1re = x + y + z + w;
                    c2re = x - y - z + w;

                    x = B0[l, m];
                    y = B0[l, m + M2];
                    z = B0[l, M2 - m].Conjugate;
                    w = B0[l, Maths.Mod(M - m, M)].Conjugate;

                    c1im = x + y - z - w;
                    c2im = x - y + z - w;

                    A1[l, m] = 1.0 / (2) * (c1re + Maths.I * c2re * exp[m]);
                    B1[l, m] = 1.0 / (2 * Maths.I) * (c1im + Maths.I * c2im * exp[m]);
                }
            }

            A1 = FFT.Backward(Matrice.Conjugate(A1));
            B1 = FFT.Backward(Matrice.Conjugate(B1));

            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    i = l * M + 2 * k;
                    j = l * M + 2 * k + 1;

                    x = A1[l, k];
                    y = B1[l, k];

                    output[i] = x.Real + Maths.I * y.Real;
                    output[j] = x.Imag + Maths.I * y.Imag;
                }
            }

            return output;
        }
        /// <summary>
        /// Fast backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
        public static Complex32[] IWHT(Complex32[] input, float[] g0, int M)
        {
            // The function implements a fast Weil-Heisenberg direct transformation algorithm,
            // stated in the following articles:
            // A. Vahlin, "EFFICIENT ALGORITHMS FOR MODULATION AND DEMODULATION IN OFDM-SYSTEMS" [1].
            // V.M. Asiryan, V.P. Volchkov, "EFFECTIVE IMPLEMENTATION OF THE DIRECT TRANSFORMATION OF WEIL-HEISENBERG" [2].
            // The algorithm is computationally efficient for large M.

            int N = input.Length, L = N / M, M2 = M / 2, M4 = M2 / 2;
            Complex32[] output = new Complex32[N];
            Complex32[,] A1 = new Complex32[L, M];
            Complex32[,] B1 = new Complex32[L, M];
            Complex32[] exp = FastWeylHeisenbergTransform.GetRotation(M);
            Complex32 s;
            int n, k, l;

            for (k = 0; k < M; k++)
            {
                for (l = 0; l < L; l++)
                {
                    s = input[k + l * M];

                    A1[l, k] = s.Real;
                    B1[l, k] = s.Imag;
                }
            }

            Complex32[,] Za = new Complex32[L, M2];
            Complex32[,] Zb = new Complex32[L, M2];

            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    Za[l, k] = A1[l, k * 2] + Maths.I * A1[l, k * 2 + 1];
                    Zb[l, k] = B1[l, k * 2] + Maths.I * B1[l, k * 2 + 1];
                }
            }

            Za = Matrice.Conjugate(FFT.Backward(Za));
            Zb = Matrice.Conjugate(FFT.Backward(Zb));

            Complex32 a0, a1, b0, b1;
            Complex32 x, y, u, v;

            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    a0 = Za[l, k]; a1 = Za[l, Maths.Mod(M - k, M2)].Conjugate;

                    x = 1.0 / (2) * (a0 + a1);
                    y = 1.0 / (2 * Maths.I) * (a0 - a1);
                    y *= exp[k];

                    A1[l, k] = x + y;
                    A1[l, k + M2] = x - y;

                    b0 = Zb[l, k]; b1 = Zb[l, Maths.Mod(M - k, M2)].Conjugate;

                    u = 1.0 / (2) * (b0 + b1);
                    v = 1.0 / (2 * Maths.I) * (b0 - b1);
                    v *= exp[k];

                    B1[l, k] = u + v;
                    B1[l, k + M2] = u - v;
                }
            }

            for (l = 0; l < L; l++)
            {
                for (n = 0; n < N; n++)
                {
                    output[n] += A1[l, Maths.Mod(n - M4, M)] * g0[Maths.Mod(n - l * M, N)] - Maths.I
                               * B1[l, Maths.Mod(n - M4, M)] * g0[Maths.Mod(n - l * M + M2, N)];
                }
            }

            return output;
        }
        /// <summary>
        /// Returns an array of phase rotations.
        /// </summary>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
        private static Complex32[] GetRotation(int M)
        {
            int M2 = M / 2;
            Complex32[] phase = new Complex32[M2];
            for (int k = 0; k < M2; k++)
            {
                phase[k] = Maths.Exp(Maths.I * 2 * Math.PI / M * k);
            }
            return phase;
        }
        #endregion
    }
}
