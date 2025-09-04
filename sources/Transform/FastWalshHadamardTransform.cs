using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Walsh-Hadamard transform.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.mathworks.com/matlabcentral/fileexchange/6879-fast-walsh-hadamard-transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastWalshHadamardTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        private bool normalized;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Walsh-Hadamard transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastWalshHadamardTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.normalized = normalized; this.direction = direction;
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
        #endregion

        #region Fast Walsh-Hadamard Transform
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            int N = A.Length;
            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[] B = (float[])A.Clone();
            FWHT(B);

            if (normalized)
            {
                B = Matrix.Div(B, Maths.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length;
            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[] A = (float[])B.Clone();
            FWHT(A);

            if (normalized)
            {
                A = Matrix.Div(A, Maths.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            float[,] B = (float[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    float[] row = new float[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    FWHT(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    float[] col = new float[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    FWHT(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrix.Div(B, Maths.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    float[] col = new float[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    FWHT(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrix.Div(B, Maths.Sqrt(N));
                }
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    float[] row = new float[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    FWHT(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    B = Matrix.Div(B, Maths.Sqrt(M));
                }
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            float[,] A = (float[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    float[] col = new float[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    FWHT(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    float[] row = new float[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    FWHT(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrix.Div(A, Maths.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    float[] col = new float[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    FWHT(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrix.Div(A, Maths.Sqrt(N));
                }
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    float[] row = new float[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    FWHT(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    A = Matrix.Div(A, Maths.Sqrt(M));
                }
            }

            return A;
        }
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            Complex32[] B = (Complex32[])A.Clone();
            FWHT(B);

            if (normalized)
            {
                B = Matrix.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            Complex32[] A = (Complex32[])B.Clone();
            FWHT(A);

            if (normalized)
            {
                A = Matrix.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            Complex32[,] B = (Complex32[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

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

                    FWHT(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    FWHT(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrix.Div(B, Math.Sqrt(N * M));
                }
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

                    FWHT(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrix.Div(B, Math.Sqrt(N));
                }
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

                    FWHT(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    B = Matrix.Div(B, Math.Sqrt(M));
                }
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
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
                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    FWHT(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    FWHT(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrix.Div(A, Math.Sqrt(N * M));
                }
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
                    FWHT(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrix.Div(A, Math.Sqrt(N));
                }
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
                    FWHT(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    A = Matrix.Div(A, Math.Sqrt(M));
                }
            }

            return A;
        }
        #endregion

        #region Private data
        /// <summary>
        /// Fast Walsh-Hadamard transform.
        /// </summary>
        /// <param name="data">Array</param>
        private void FWHT(float[] data)
        {
            int N = data.Length;

            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int log2N = (int)Maths.Log2(N);
            float x_even, x_odd;

            int k0 = N, k1 = 1, k2 = N / 2;
            int x, y, z, i, j, l;

            for (x = 0; x < log2N; x++)
            {
                l = 0;

                for (y = 0; y < k1; y++, l += k0)
                {
                    for (z = 0; z < k2; z++)
                    {
                        i = z + l; j = i + k2;

                        x_even = data[i];
                        x_odd = data[j];

                        data[i] = x_even + x_odd;
                        data[j] = x_even - x_odd;
                    }
                }

                k0 /= 2; k1 *= 2; k2 /= 2;
            }
        }
        /// <summary>
        /// Fast Walsh-Hadamard transform.
        /// </summary>
        /// <param name="data">Array</param>
        private void FWHT(Complex32[] data)
        {
            int N = data.Length;

            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int log2N = (int)Maths.Log2(N);
            Complex32 x_even, x_odd;

            int k0 = N, k1 = 1, k2 = N / 2;
            int x, y, z, i, j, l;

            for (x = 0; x < log2N; x++)
            {
                l = 0;

                for (y = 0; y < k1; y++, l += k0)
                {
                    for (z = 0; z < k2; z++)
                    {
                        i = z + l; j = i + k2;

                        x_even = data[i];
                        x_odd = data[j];

                        data[i] = x_even + x_odd;
                        data[j] = x_even - x_odd;
                    }
                }

                k0 /= 2; k1 *= 2; k2 /= 2;
            }
        }
        #endregion
    }
}
