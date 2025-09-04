using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Hilbert transform.
    /// <remarks>
    /// NOT RECOMMENDED.
    /// 
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
        private readonly FastFourierTransform FFT;
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
        public Complex32[] Forward(Complex32[] A)
        {
            var F = FFT.Forward(A);
            HilbertTransform.ApplyHilbertOperatorMaskInplace(F);
            return FFT.Backward(F);
        }
        /// <summary>
        /// Backward Hilbert transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            var Hx = Forward(B);
            for (int i = 0; i < Hx.Length; i++) Hx[i] = -Hx[i];
            return Hx;
        }
        /// <summary>
        /// Forward Hilbert transform.
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

                    row = Forward(row);

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
        /// Backward Hilbert transform.
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
        /// Forward Hilbert transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Hilbert transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Hilbert transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Hilbert transform.
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
