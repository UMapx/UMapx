using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Laplace transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Laplace_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastLaplaceTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FastFourierTransform FFT;
        /// <summary>
        /// Damping factor.
        /// </summary>
        private float sigma;
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Laplace transform.
        /// </summary>
        /// <param name="sigma">Damping factor (0, 1)</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastLaplaceTransform(float sigma = 0.0005f, bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(normalized, Direction.Vertical);
            this.Sigma = sigma; 
            this.Direction = direction;
        }
        /// <summary>
        /// Gets or sets the damping factor (0, 1).
        /// <remarks>
        /// If σ = 0, then the Laplace transform takes the form of a Fourier transform.
        /// </remarks>
        /// </summary>
        public float Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.sigma = value;
            }
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

        #region Laplace Transform
        /// <summary>
        /// Forward Laplace transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public ComplexF[] Forward(ComplexF[] A)
        {
            // Fourier transform:
            ComplexF[] B = FFT.Forward(A);

            // Fourier to Laplace transform:
            ApplyLaplaceOperator(B, sigma);
            return B;
        }
        /// <summary>
        /// Backward Laplace transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public ComplexF[] Backward(ComplexF[] B)
        {
            // Laplace to Fourier transform:
            ComplexF[] A = (ComplexF[])B.Clone();
            ApplyLaplaceInverseOperator(A, sigma);

            // Fourier transform:
            return FFT.Backward(A);
        }
        /// <summary>
        /// Forward Laplace transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public ComplexF[,] Forward(ComplexF[,] A)
        {
            ComplexF[,] B = (ComplexF[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    ComplexF[] row = new ComplexF[M];
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
                    ComplexF[] col = new ComplexF[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j].Conjugate;
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
                    ComplexF[] col = new ComplexF[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j].Conjugate;
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
                    ComplexF[] row = new ComplexF[M];
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
        /// Backward Laplace transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public ComplexF[,] Backward(ComplexF[,] B)
        {
            ComplexF[,] A = (ComplexF[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    ComplexF[] col = new ComplexF[N];
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
                    ComplexF[] row = new ComplexF[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j].Conjugate;
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
                    ComplexF[] col = new ComplexF[N];
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
                    ComplexF[] row = new ComplexF[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j].Conjugate;
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
        /// Forward Laplace transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Laplace transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Laplace transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Laplace transform.
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
        /// Applies Laplace transform operator.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="sigma">Sigma</param>
        private static void ApplyLaplaceOperator(ComplexF[] v, float sigma = 0.003f)
        {
            // forward laplacian transform
            // for signal:
            int N = v.Length;

            for (int i = 0; i < N; i++)
            {
                v[i] = Math.Exp(-sigma * i) * v[i];
            }
        }
        /// <summary>
        /// Applies inverse Laplace transform operator.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="sigma">Sigma</param>
        private static void ApplyLaplaceInverseOperator(ComplexF[] v, float sigma = 0.003f)
        {
            // inverse laplacian transform
            // for signal:
            int N = v.Length;

            for (int i = 0; i < N; i++)
            {
                v[i] = v[i] / Math.Exp(-sigma * i);
            }
        }
        #endregion
    }
}
