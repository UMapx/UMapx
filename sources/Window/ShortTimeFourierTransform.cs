using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines short-time Fourier transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Short-time_Fourier_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ShortTimeFourierTransform : IWindowTransform, ITransform
    {
        #region Private data
        private readonly FourierTransform DFT;
        private IWindow window;
        private Direction direction;
        private readonly float[] coefs;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes short-time Fourier transform.
        /// </summary>
        /// <param name="function">Windows function</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public ShortTimeFourierTransform(IWindow function, bool normalized = true, Direction direction = Direction.Vertical)
        {
            // fourier transform initialization:
            this.DFT = new FourierTransform(normalized, direction);
            Direction = direction;
            Window = function;

            // sampling window function:
            this.coefs = function.GetWindow().Add(1e-16f);
        }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.DFT.Normalized;
            }
            set
            {
                this.DFT.Normalized = value;
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

        #region Short-time Fourier transform
        /// <summary>
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            // params
            int N = A.Length, i, j;
            Complex32[] B = new Complex32[N];
            int frame = coefs.Length;

            // Short-Time Fourier Transform
            for (i = 0; i < N; i += frame)
            {
                Complex32[] data = new Complex32[frame];

                for (j = 0; j < frame; j++)
                    data[j] = A[i + j] * coefs[Maths.Mod(i - frame / 2, frame)];

                data = DFT.Forward(data);

                for (j = 0; j < frame; j++)
                    B[i + j] = data[j];
            }

            return B;
        }
        /// <summary>
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length, i, j;
            Complex32[] A = new Complex32[N];
            int frame = coefs.Length;

            for (i = 0; i < N; i += frame)
            {
                Complex32[] data = new Complex32[frame];

                for (j = 0; j < frame; j++)
                {
                    data[j] = B[i + j];
                }

                data = DFT.Backward(data);

                for (j = 0; j < frame; j++)
                {
                    A[i + j] = data[j] / coefs[Maths.Mod(i - frame / 2, frame)];
                }
            }

            return A;
        }
        /// <summary>
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            // Fourier transform:
            Complex32[,] B = (Complex32[,])A.Clone();
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            Complex32[,] A = (Complex32[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
}
