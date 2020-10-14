using System;
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
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            if (!Maths.IsPower(N, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            Complex[] B = (Complex[])A.Clone();
            fft(B);

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
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            if (!Maths.IsPower(N, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            Complex[] A = (Complex[])B.Clone();
            ifft(A);

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
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = (Complex[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                    throw new Exception("Dimension of the signal must be a power of 2");

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    ifft(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    fft(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {
                if (!Maths.IsPower(N, 2))
                    throw new Exception("Dimension of the signal must be a power of 2");

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    fft(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N));
                }
            }
            else
            {
                if (!Maths.IsPower(M, 2))
                    throw new Exception("Dimension of the signal must be a power of 2");

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    ifft(row);

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
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                    throw new Exception("Dimension of the signal must be a power of 2");

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    ifft(col);

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
                    fft(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {
                if (!Maths.IsPower(N, 2))
                    throw new Exception("Dimension of the signal must be a power of 2");

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    ifft(col);

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
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    fft(row);

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

        #region Private data
        private const int minLength = 2;
        private const int maxLength = 16384;
        private const int minBits = 1;
        private const int maxBits = 14;
        private static int[][] reversedBits = new int[maxBits][];
        private static Complex[,][] complexRotation = new Complex[maxBits, 2][];
        #endregion

        #region Private voids
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="data">Array</param>
        private static void fft(Complex[] data)
        {
            int n = data.Length;
            int m = Log2(n);

            // reorder data first
            FastFourierTransform.ReorderData(data);

            // compute FFT
            int tn = 1, tm, k, i, even, odd;
            Complex[] rotation;
            Complex t, ce, co;
            double tr, ti;

            for (k = 1; k <= m; k++)
            {
                rotation = FastFourierTransform.ForwardComplexRotation(k);
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
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="data">Array</param>
        private static void ifft(Complex[] data)
        {
            int n = data.Length;
            int m = Log2(n);

            // reorder data first
            FastFourierTransform.ReorderData(data);

            // compute FFT
            int tn = 1, tm, k, i, even, odd;
            Complex[] rotation;
            Complex t, ce, co;
            double tr, ti;

            for (k = 1; k <= m; k++)
            {
                rotation = FastFourierTransform.BackwardComplexRotation(k);
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
        private static Complex[] ForwardComplexRotation(int numberOfBits)
        {
            int directionIndex = 0;

            // check if the array is already calculated
            if (complexRotation[numberOfBits - 1, directionIndex] == null)
            {
                int n = 1 << (numberOfBits - 1), i;
                double uR = 1.0;
                double uI = 0.0;
                double angle = -System.Math.PI / n;
                double wR = System.Math.Cos(angle);
                double wI = System.Math.Sin(angle);
                double t;
                Complex[] rotation = new Complex[n];

                for (i = 0; i < n; i++)
                {
                    rotation[i] = new Complex(uR, uI);
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
        private static Complex[] BackwardComplexRotation(int numberOfBits)
        {
            int directionIndex = 1;

            // check if the array is already calculated
            if (complexRotation[numberOfBits - 1, directionIndex] == null)
            {
                int n = 1 << (numberOfBits - 1), i;
                double uR = 1.0;
                double uI = 0.0;
                double angle = System.Math.PI / n;
                double wR = System.Math.Cos(angle);
                double wI = System.Math.Sin(angle);
                double t;
                Complex[] rotation = new Complex[n];

                for (i = 0; i < n; i++)
                {
                    rotation[i] = new Complex(uR, uI);
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
        private static void ReorderData(Complex[] data)
        {
            int length = data.Length;
            int[] rBits = GetReversedBits(Log2(length));
            Complex t;
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
            return (int)(Math.Log10(x) / 0.30102999566398);
        }
        #endregion
    }
}
