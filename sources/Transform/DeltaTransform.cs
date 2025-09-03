using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the delta transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Delta_encoding
    /// </remarks>
    /// </summary>
    [Serializable]
    public class DeltaTransform : ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the delta transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public DeltaTransform(Direction direction = Direction.Vertical)
        {
            Direction = direction;
        }
        /// <summary>
        /// Gets or sets the processing direction.
        /// </summary>
        public Direction Direction { get; set; }
        #endregion

        #region Delta transform
        /// <summary>
        /// Forward delta transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            float[] B = (float[])A.Clone();
            int length = B.Length;
            float last = 0;

            for (int i = 0; i < length; i++)
            {
                var current = B[i];
                B[i] = current - last;
                last = current;
            }

            return B;
        }
        /// <summary>
        /// Backward delta transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            float[] A = (float[])B.Clone();
            int length = A.Length;
            float last = 0;

            for (int i = 0; i < length; i++)
            {
                var delta = A[i];
                A[i] = delta + last;
                last = A[i];
            }

            return A;
        }
        /// <summary>
        /// Forward delta transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            float[,] B = (float[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (Direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    float[] row = new float[M];
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
                    float[] col = new float[N];
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
            else if (Direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    float[] col = new float[N];
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
                    float[] row = new float[M];
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
        /// Backward delta transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            float[,] A = (float[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (Direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    float[] col = new float[N];
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
                    float[] row = new float[M];
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
                Parallel.For(0, M, j =>
                {
                    float[] col = new float[N];
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
                    float[] row = new float[M];
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
        /// Forward delta transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public ComplexF[] Forward(ComplexF[] A)
        {
            ComplexF[] B = (ComplexF[])A.Clone();
            int length = B.Length;
            ComplexF last = 0;

            for (int i = 0; i < length; i++)
            {
                var current = B[i];
                B[i] = current - last;
                last = current;
            }

            return B;
        }
        /// <summary>
        /// Backward delta transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public ComplexF[] Backward(ComplexF[] B)
        {
            ComplexF[] A = (ComplexF[])B.Clone();
            int length = A.Length;
            ComplexF last = 0;

            for (int i = 0; i < length; i++)
            {
                var delta = A[i];
                A[i] = delta + last;
                last = A[i];
            }

            return A;
        }
        /// <summary>
        /// Forward delta transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public ComplexF[,] Forward(ComplexF[,] A)
        {
            ComplexF[,] B = (ComplexF[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (Direction == Direction.Both)
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
            else if (Direction == Direction.Vertical)
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
        /// Backward delta transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public ComplexF[,] Backward(ComplexF[,] B)
        {
            ComplexF[,] A = (ComplexF[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (Direction == Direction.Both)
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
                }
                );
            }
            else if (Direction == Direction.Vertical)
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
        #endregion
    }
}
