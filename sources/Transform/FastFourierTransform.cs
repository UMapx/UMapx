using System;
using System.Numerics;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Fourier transform using the Cooley-Tukey algorithm.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastFourierTransform : ITransform
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
        /// Initializes the fast Fourier transform using the Cooley-Tukey algorithm.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastFourierTransform(bool normalized = true, Direction direction = Direction.Vertical)
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

        #region Fast Fourier Transform
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            Complex32[] B = (Complex32[])A.Clone();

            if (!Maths.IsPower(N, 2))
            {
                BluesteinFFT(B, false);
            }
            else
            {
                CooleyTukeyFFT(B, false);
            }

            if (normalized == true)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            Complex32[] A = (Complex32[])B.Clone();

            if (!Maths.IsPower(N, 2))
            {
                BluesteinFFT(B, true);
            }
            else
            {
                CooleyTukeyFFT(A, true);
            }

            if (normalized == true)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Fourier transform.
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

                    if (!Maths.IsPower(M, 2))
                    {
                        BluesteinFFT(row, true);
                    }
                    else
                    {
                        CooleyTukeyFFT(row, true);
                    }

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                Parallel.For(0, M, j =>
                {
                    Complex32[] col = new Complex32[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    if (!Maths.IsPower(N, 2))
                    {
                        BluesteinFFT(col, false);
                    }
                    else
                    {
                        CooleyTukeyFFT(col, false);
                    }

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N * M));
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

                    if (!Maths.IsPower(N, 2))
                    {
                        BluesteinFFT(col, false);
                    }
                    else
                    {
                        CooleyTukeyFFT(col, false);
                    }

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N));
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

                    if (!Maths.IsPower(M, 2))
                    {
                        BluesteinFFT(row, true);
                    }
                    else
                    {
                        CooleyTukeyFFT(row, true);
                    }

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(M));
                }
            }

            return B;
        }
        /// <summary>
        /// Backward Fourier transform.
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

                    if (!Maths.IsPower(N, 2))
                    {
                        BluesteinFFT(col, true);
                    }
                    else
                    {
                        CooleyTukeyFFT(col, true);
                    }

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    if (!Maths.IsPower(M, 2))
                    {
                        BluesteinFFT(row, false);
                    }
                    else
                    {
                        CooleyTukeyFFT(row, false);
                    }

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(N * M));
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

                    if (!Maths.IsPower(N, 2))
                    {
                        BluesteinFFT(col, true);
                    }
                    else
                    {
                        CooleyTukeyFFT(col, true);
                    }

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(N));
                }
            }
            else
            {
                if (!Maths.IsPower(M, 2))
                    throw new Exception("Dimension of the signal must be a power of 2");

                Parallel.For(0, N, i =>
                {
                    Complex32[] row = new Complex32[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    if (!Maths.IsPower(M, 2))
                    {
                        BluesteinFFT(row, false);
                    }
                    else
                    {
                        CooleyTukeyFFT(row, false);
                    }

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(M));
                }
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

        #region Private data
        private const int minBits = 1;
        private const int maxBits = 14;
        private static int[][] reversedBits = new int[maxBits][];
        private static Complex32[,][] complexRotation = new Complex32[maxBits, 2][];
        #endregion

        #region Private voids
        /// <summary>
        /// Forward Fourier transform (Bluestein FFT).
        /// </summary>
        /// <param name="input"></param>
        /// <param name="inverse"></param>
        public static void BluesteinFFT(Complex32[] input, bool inverse)
        {
            int N = input.Length;

            // Находим ближайшую степень двойки
            int M = 1;
            while (M < 2 * N - 1)
                M <<= 1;

            var a = new Complex32[M];
            var b = new Complex32[M];
            var c = new Complex32[M];

            float sign = inverse ? 1f : -1f;
            float norm = 1f / M;

            // Chirp-модуляция
            for (int n = 0; n < N; n++)
            {
                float angle = Maths.Pi * n * n / N;
                Complex32 w = Maths.Exp(sign * Maths.I * angle);
                a[n] = input[n] * w;
                b[n] = Maths.Exp(-sign * Maths.I * angle);
            }

            // Остальное обнулено
            for (int i = N; i < M; i++)
            {
                a[i] = Complex32.Zero;
                b[i] = Complex32.Zero;
            }

            // Зеркальное дополнение ядра
            for (int i = 1; i < N; i++)
                b[M - i] = b[i];

            // Выполняем свёртку через FFT
            CooleyTukeyFFT(a, false);
            CooleyTukeyFFT(b, false);
            for (int i = 0; i < M; i++)
                c[i] = a[i] * b[i];
            CooleyTukeyFFT(c, true);

            // Обратная модуляция + нормализация
            for (int n = 0; n < N; n++)
            {
                float angle = Maths.Pi * n * n / N;
                Complex32 w = Maths.Exp(sign * Maths.I * angle);
                input[n] = c[n] * w;
                input[n] *= norm; // Нормализация для согласованности с матрицей
            }
        }
        /// <summary>
        /// Forward Fourier transform (Cooley-Tukey FFT).
        /// </summary>
        /// <param name="data">Array</param>
        /// <param name="inverse">Inverse or not</param>
        private static void CooleyTukeyFFT(Complex32[] data, bool inverse)
        {
            int n = data.Length;
            int m = Log2(n);

            // reorder data first
            FastFourierTransform.ReorderData(data);

            // compute FFT
            int tn = 1, tm, k, i, even, odd;
            Complex32[] rotation;
            Complex32 t, ce, co;
            float tr, ti;

            for (k = 1; k <= m; k++)
            {
                rotation = inverse ? FastFourierTransform.BackwardComplexRotation(k) : FastFourierTransform.ForwardComplexRotation(k);
                tm = tn; tn <<= 1;

                for (i = 0; i < tm; i++)
                {
                    t = rotation[i];

                    for (even = i; even < n; even += tn)
                    {
                        odd = even + tm;
                        ce = data[even];
                        co = data[odd];

                        tr = co.Real * t.Real - co.Imag * t.Imag;
                        ti = co.Real * t.Imag + co.Imag * t.Real;

                        data[even].Real += tr;
                        data[even].Imag += ti;

                        data[odd].Real = ce.Real - tr;
                        data[odd].Imag = ce.Imag - ti;
                    }
                }
            }
        }
        /// <summary>
        /// Gets an array with pointers to data members that must be replaced before the FFT.
        /// </summary>
        /// <param name="numberOfBits">Number of bits</param>
        /// <returns>Array</returns>
        private static int[] GetReversedBits(int numberOfBits)
        {
            if ((numberOfBits < minBits) || (numberOfBits > maxBits))
                throw new ArgumentOutOfRangeException();

            // check if the array is already calculated
            if (reversedBits[numberOfBits - 1] == null)
            {
                int n = Pow2(numberOfBits);
                int[] rBits = new int[n];
                int i, j, oldBits, newBits;

                // calculate the array
                for (i = 0; i < n; i++)
                {
                    oldBits = i;
                    newBits = 0;

                    for (j = 0; j < numberOfBits; j++)
                    {
                        newBits = (newBits << 1) | (oldBits & 1);
                        oldBits = (oldBits >> 1);
                    }
                    rBits[i] = newBits;
                }
                reversedBits[numberOfBits - 1] = rBits;
            }
            return reversedBits[numberOfBits - 1];
        }
        /// <summary>
        /// Gets the forward rotation of a complex number.
        /// </summary>
        /// <param name="numberOfBits">Number of bits</param>
        /// <returns>Array</returns>
        private static Complex32[] ForwardComplexRotation(int numberOfBits)
        {
            int directionIndex = 0;

            // check if the array is already calculated
            if (complexRotation[numberOfBits - 1, directionIndex] == null)
            {
                int n = 1 << (numberOfBits - 1), i;
                float uR = 1.0f;
                float uI = 0.0f;
                float angle = -Maths.Pi / n;
                float wR = Maths.Cos(angle);
                float wI = Maths.Sin(angle);
                float t;
                Complex32[] rotation = new Complex32[n];

                for (i = 0; i < n; i++)
                {
                    rotation[i] = new Complex32(uR, uI);
                    t = uR * wI + uI * wR;
                    uR = uR * wR - uI * wI;
                    uI = t;
                }

                complexRotation[numberOfBits - 1, directionIndex] = rotation;
            }
            return complexRotation[numberOfBits - 1, directionIndex];
        }
        /// <summary>
        /// Gets the backward rotation of a complex number.
        /// </summary>
        /// <param name="numberOfBits">Number of bits</param>
        /// <returns>Array</returns>
        private static Complex32[] BackwardComplexRotation(int numberOfBits)
        {
            int directionIndex = 1;

            // check if the array is already calculated
            if (complexRotation[numberOfBits - 1, directionIndex] == null)
            {
                int n = 1 << (numberOfBits - 1), i;
                float uR = 1.0f;
                float uI = 0.0f;
                float angle = Maths.Pi / n;
                float wR = Maths.Cos(angle);
                float wI = Maths.Sin(angle);
                float t;
                Complex32[] rotation = new Complex32[n];

                for (i = 0; i < n; i++)
                {
                    rotation[i] = new Complex32(uR, uI);
                    t = uR * wI + uI * wR;
                    uR = uR * wR - uI * wI;
                    uI = t;
                }

                complexRotation[numberOfBits - 1, directionIndex] = rotation;
            }
            return complexRotation[numberOfBits - 1, directionIndex];
        }
        /// <summary>
        /// Reorders data to use FFT.
        /// </summary>
        /// <param name="data">Array</param>
        private static void ReorderData(Complex32[] data)
        {
            int length = data.Length;
            int[] rBits = GetReversedBits(Log2(length));
            Complex32 t;
            int i, s;

            for (i = 0; i < length; i++)
            {
                s = rBits[i];

                if (s > i)
                {
                    t = data[i];
                    data[i] = data[s];
                    data[s] = t;
                }
            }
        }
        /// <summary>
        /// Computes power of 2.
        /// </summary>
        /// <param name="power">Power</param>
        /// <returns>Integer number</returns>
        private static int Pow2(int power)
        {
            return ((power >= 0) && (power <= 30)) ? (1 << power) : 0;
        }
        /// <summary>
        /// Calculates the base 2 logarithm.
        /// </summary>
        /// <param name="x">Integer number</param>
        /// <returns>Integer number</returns>
        private static int Log2(int x)
        {
            return (int)Math.Log(x, 2);
        }
        #endregion
    }
}
