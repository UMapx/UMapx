using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Hartley transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastHartleyTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FastFourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Hartley transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastHartleyTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(normalized, direction);
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

        #region Fast Hartley Transform
        /// <summary>
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            Complex32[] B = Matrice.ToComplex(A);
            B = FFT.Forward(B);

            int length = A.Length, i;
            float[] Hk = new float[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = B[i].Real - B[i].Imag;
            }

            return Hk;
        }
        /// <summary>
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            Complex32[] A = Matrice.ToComplex(B);
            A = FFT.Backward(A);

            int length = B.Length, i;
            float[] Hk = new float[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = A[i].Real + A[i].Imag;
            }

            return Hk;
        }
        /// <summary>
        /// Forward Hartley transform.
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
        /// Backward Hartley transform.
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
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int n = A.Length;
            var re = new float[n];
            var im = new float[n];

            for (int i = 0; i < n; i++)
            {
                re[i] = A[i].Real;
                im[i] = A[i].Imag;
            }

            // reuse correct real FHT you showed above
            var Hr = Forward(re);
            var Hi = Forward(im);

            var H = new Complex32[n];
            for (int i = 0; i < n; i++)
                H[i] = new Complex32(Hr[i], Hi[i]);

            return H;
        }
        /// <summary>
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int n = B.Length;
            var re = new float[n];
            var im = new float[n];

            for (int i = 0; i < n; i++)
            {
                re[i] = B[i].Real;
                im[i] = B[i].Imag;
            }

            var xr = Backward(re);
            var xi = Backward(im);

            var X = new Complex32[n];
            for (int i = 0; i < n; i++)
                X[i] = new Complex32(xr[i], xi[i]);

            return X;
        }
        /// <summary>
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            // split into two real planes
            var R = new float[N, M];
            var I = new float[N, M];

            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                {
                    R[i, j] = A[i, j].Real;
                    I[i, j] = A[i, j].Imag;
                }

            // reuse your correct real 2D FHT
            var HR = Forward(R);
            var HI = Forward(I);

            // combine back into complex result
            var H = new Complex32[N, M];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                    H[i, j] = new Complex32(HR[i, j], HI[i, j]);

            return H;
        }
        /// <summary>
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            var R = new float[N, M];
            var I = new float[N, M];

            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                {
                    R[i, j] = B[i, j].Real;
                    I[i, j] = B[i, j].Imag;
                }

            var XR = Backward(R);
            var XI = Backward(I);

            var X = new Complex32[N, M];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                    X[i, j] = new Complex32(XR[i, j], XI[i, j]);

            return X;
        }

        #endregion
    }
}
