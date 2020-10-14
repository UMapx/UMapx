using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast sine transform.
    /// <remarks>
    /// More information can be found on the website:
    /// http://sernam.ru/book_prett1.php?id=91
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastSineTransform : ITransform
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
        /// Defines the fast sine transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public FastSineTransform(Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(true, Direction.Both);
            this.direction = direction;
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

        #region Fast Sine Transform
        /// <summary>
        /// Forward cosine transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length, N2 = N * 2, i, k;
            Complex[] B = new Complex[N2];

            for (i = 0; i < N; i++)
            {
                B[i] = A[i];
            }

            B = FFT.Forward(B);

            double[] C = new double[N];
            Complex c = -Maths.I * Maths.Pi / N;

            for (k = 0; k < N; k++)
            {
                C[k] = 2.0 * (B[k] * Maths.Exp(c * k)).Imag;
            }

            C[0] = A[N - 1] / Math.Sqrt(2);

            return C;
        }
        /// <summary>
        /// Backward cosine transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length, N2 = N * 2, i, k;
            Complex[] A = new Complex[N2];
            double Bk, temp, c = Maths.Pi / N;

            for (k = 0; k < N; k++)
            {
                Bk = B[k];
                temp = k * c;
                A[k] = new Complex(Bk * Math.Cos(temp), Bk * Math.Sin(temp));
            }

            A = FFT.Backward(A);
            double[] C = new double[N];

            for (i = 0; i < N; i++)
            {
                C[i] = -A[i].Imag * 2.0;
            }

            C[N - 1] = B[0] * Math.Sqrt(2);

            return C;
        }
        /// <summary>
        /// Forward cosine transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            double[,] B = (double[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
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
                }
                );

                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
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
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
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
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
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
        /// Backward cosine transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            double[,] A = (double[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
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
                    double[] row = new double[M];
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
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
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
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
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
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
}
