using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Hilbert transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hilbert_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastHilbertTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private FastFourierTransform FFT;
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Hilbert transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastHilbertTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(normalized, Direction.Both);
            this.direction = direction;
        }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
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

        #region Fast Hilbert Transform
        /// <summary>
        /// Forward Hilbert transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            Complex[] F = FFT.Forward(A);
            HilbertTransform.hilbertf(F, N);
            F = FFT.Backward(F);
            return HilbertTransform.hilbertb(A, F, N);
        }
        /// <summary>
        /// Backward Hilbert transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length, i;
            Complex[] A = new Complex[N];

            for (i = 0; i < N; i++)
            {
                A[i] = new Complex(B[i].Real, 0);
            }

            return A;
        }
        /// <summary>
        /// Forward Hilbert transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = (Complex[,])A.Clone();
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            if (Direction == Direction.Both)
            {
                // 2-dimension horizontal Hilbert transform:
                FFT.Direction = Direction.Horizontal;
                B = FFT.Forward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    HilbertTransform.hilbertf(row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    Complex[] num = new Complex[M];

                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                        num[j] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = num[j];
                    }
                });


                // 2-dimension vertical Hilbert transform:
                FFT.Direction = Direction.Vertical;
                B = FFT.Forward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    HilbertTransform.hilbertf(col, N);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    Complex[] num = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                        num[i] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, col, N);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = num[i];
                    }
                });

                FFT.Direction = Direction.Both;
            }
            else if (Direction == Direction.Vertical)
            {
                // 2-dimension vertical Hilbert transform:
                B = FFT.Forward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    HilbertTransform.hilbertf(col, N);

                    for (i = 0; i < M; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    Complex[] num = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                        num[i] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, col, N);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = num[i];
                    }
                });
            }
            else
            {
                // 2-dimension horizontal Hilbert transform:
                B = FFT.Forward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    HilbertTransform.hilbertf(row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    Complex[] num = new Complex[M];

                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                        num[j] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = num[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Backward Hilbert transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (Direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
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
                }
                );

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
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
                }
                );
            }
            else if (Direction == Direction.Vertical)
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
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
                }
                );
            }
            else
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
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

            return A;
        }
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
}
