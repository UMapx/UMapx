using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Hilbert transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hilbert_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class HilbertTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private FourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Hilbert transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public HilbertTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FourierTransform(normalized, direction);
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
                return this.FFT.Direction;
            }
            set
            {
                this.FFT.Direction = value;
            }
        }
        #endregion

        #region Hilbert Transform
        /// <summary>
        /// Forward Hilbert transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            Complex32[] F = FFT.Forward(A);
            HilbertTransform.hilbertf(F, N);
            F = FFT.Backward(F);
            return HilbertTransform.hilbertb(A, F, N);
        }
        /// <summary>
        /// Backward Hilbert transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length, i;
            Complex32[] A = new Complex32[N];

            for (i = 0; i < N; i++)
            {
                A[i] = new Complex32(B[i].Real, 0);
            }

            return A;
        }
        /// <summary>
        /// Forward Hilbert transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            Complex32[,] B = (Complex32[,])A.Clone();
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            if (Direction == Direction.Both)
            {
                // 2-dimension horizontal Hilbert transform:
                FFT.Direction = Direction.Horizontal;
                B = FFT.Forward(B);

                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
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
                    Complex32[] row = new Complex32[M];
                    Complex32[] num = new Complex32[M];

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
                    Complex32[] col = new Complex32[N];
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
                    Complex32[] col = new Complex32[N];
                    Complex32[] num = new Complex32[N];
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
                    Complex32[] col = new Complex32[N];
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
                    Complex32[] col = new Complex32[N];
                    Complex32[] num = new Complex32[N];
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
                    Complex32[] row = new Complex32[M];
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
                    Complex32[] row = new Complex32[M];
                    Complex32[] num = new Complex32[M];

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
        public Complex32[,] Backward(Complex32[,] B)
        {
            Complex32[,] A = (Complex32[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (Direction == Direction.Both)
            {
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
                }
                );
            }
            else
            {
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

            return A;
        }
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Implements the rearrangement of the spectrum UMapxing to Hilbert.
        /// </summary>
        /// <param name="f">Spectrum</param>
        /// <param name="n">Length</param>
        internal static void hilbertf(Complex32[] f, int n)
        {
            int n2 = n / 2;

            for (int i = 0; i < n2; i++)
            {
                f[i] *= 2.0;
                f[i + n2] = Complex32.Zero;
            }
            return;
        }
        /// <summary>
        /// Implements the rearrangement of the spectrum UMapxing to Hilbert.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="f">Spectrum</param>
        /// <param name="n">Length</param>
        /// <returns>Array</returns>
        internal static Complex32[] hilbertb(Complex32[] a, Complex32[] f, int n)
        {
            Complex32[] B = new Complex32[n];

            for (int i = 0; i < n; i++)
            {
                B[i] = new Complex32(a[i].Real, f[i].Imag);
            }

            return B;
        }
        #endregion
    }
}
