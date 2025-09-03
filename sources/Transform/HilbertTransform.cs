using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Hilbert transform.
    /// <remarks>
    /// NOT RECOMMENDED.
    /// 
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
        private readonly FourierTransform DFT;
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Hilbert transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public HilbertTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.DFT = new FourierTransform(normalized, Direction.Both);
            this.direction = direction;
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
        #endregion

        #region Hilbert Transform
        /// <summary>
        /// Forward Hilbert transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public ComplexF[] Forward(ComplexF[] A)
        {
            var F = DFT.Forward(A);
            ApplyHilbertOperatorMaskInplace(F);
            return DFT.Backward(F);
        }
        /// <summary>
        /// Backward Hilbert transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public ComplexF[] Backward(ComplexF[] B)
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
                    ComplexF[] col = new ComplexF[N];
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
        /// Backward Hilbert transform.
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

        #region Private voids
        /// <summary>
        /// Applies Hilber spectrum operator.
        /// </summary>
        /// <param name="F">Array</param>
        internal static void ApplyHilbertOperatorMaskInplace(ComplexF[] F)
        {
            int N = F.Length;
            if (N == 0) return;

            // DC
            F[0] = ComplexF.Zero;

            if ((N & 1) == 0)
            {
                int half = N / 2;

                // Nyquist
                F[half] = ComplexF.Zero;

                // positive freqs: 1..half-1  -> -i
                for (int k = 1; k < half; k++) F[k] *= new ComplexF(0f, -1f);

                // negative freqs: half+1..N-1 -> +i
                for (int k = half + 1; k < N; k++) F[k] *= new ComplexF(0f, +1f);
            }
            else
            {
                int half = N / 2; // floor
                                  // positive: 1..half -> -i
                for (int k = 1; k <= half; k++) F[k] *= new ComplexF(0f, -1f);
                // negative: half+1..N-1 -> +i
                for (int k = half + 1; k < N; k++) F[k] *= new ComplexF(0f, +1f);
            }
        }
        #endregion
    }
}
