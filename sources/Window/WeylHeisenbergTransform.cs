using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines a group of orthogonal bases and discrete Weyl-Heisenberg transforms.
    /// <remarks>
    /// More information can be found on the website:
    /// https://elibrary.ru/item.asp?id=29767333
    /// </remarks>
    /// </summary>
    [Serializable]
    public class WeylHeisenbergTransform : IWindowTransform, ITransform
    {
        #region Private data
        /// <summary>
        /// Windows function.
        /// </summary>
        private IWindow window;
        /// <summary>
        /// Number of frequency shifts.
        /// </summary>
        private int m;
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes a group of orthogonal bases and Weyl-Heisenberg transformations.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N]</param>
        /// <param name="direction">Processing direction</param>
        public WeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical)
        {
            Window = window; M = m; Direction = direction;
        }
        /// <summary>
        /// Gets or sets number of frequency shifts [4, N].
        /// <remarks>
        /// Even number.
        /// </remarks>
        /// </summary>
        public int M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value <= 2 || !Maths.IsEven(value))
                    throw new Exception("Invalid argument value");

                this.m = value;
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
        /// <summary>
        /// Gets or sets the window function.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Weyl-Heisenberg static components
        /// <summary>
        /// Returns the complex Weyl-Heisenberg basis matrix.
        /// <remarks>
        /// Matrix dimension[N, N], where N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="N">Number of samples</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <param name="orthogonalize">Orthogonalized matrix or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] WeylHeisenberg(IWindow window, int N, int M, bool orthogonalize = true)
        {
            return WeylHeisenbergTransform.WeylHeisenberg(WeylHeisenbergTransform.GetPacket(window, N), M, orthogonalize);
        }
        /// <summary>
        /// Returns the complex Weyl-Heisenberg basis matrix.
        /// <remarks>
        /// Matrix dimension[N, N], where N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <param name="orthogonalize">Orthogonalized matrix or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] WeylHeisenberg(double[] g0, int M, bool orthogonalize = true)
        {
            if (orthogonalize)
            {
                return WeylHeisenbergTransform.WeylHeisenberg(WeylHeisenbergTransform.Zak(g0, M), M);
            }

            return WeylHeisenbergTransform.WeylHeisenberg(g0, M);
        }
        /// <summary>
        /// Returns a vector of window function values.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="length">Number of samples</param>
        /// <returns>Array</returns>
        public static double[] GetPacket(IWindow window, int length)
        {
            // exeption by length
            if (window.FrameSize > length)
                return WeylHeisenbergTransform.nSymmetry(window, length);

            // params for approximation
            double[] w = WeylHeisenbergTransform.nSymmetry(window, (int)window.FrameSize);
            int n = w.Length;
            double min = Math.Min(w[0], w[n - 1]);
            double[] g = new double[length];
            int i, j = (length - n) / 2;
            int k = Math.Min(length - 2 * j, n);
            int z = j + k;

            // do job for intervals
            for (i = 0; i < j; i++)
                g[i] = min;

            for (i = j; i < z; i++)
                g[i] = w[i - j];

            for (i = z; i < length; i++)
                g[i] = min;

            return g;
        }
        /// <summary>
        /// Returns a vector of values of a window function that satisfies the N-1 symmetry condition.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="length">Number of samples</param>
        /// <returns>Array</returns>
        private static double[] nSymmetry(IWindow window, int length)
        {
            // creaing window function
            double[] g = window.GetWindow(length + 1);
            double[] w = new double[length];

            // N-1 symmetric
            for (int i = 0; i < length; i++)
                w[i] = g[i];

            return w;
        }
        /// <summary>
        /// Returns the complex Weyl-Heisenberg basis matrix.
        /// <remarks>
        /// Matrix dimension[N, N], where N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Matrix</returns>
        private static Complex[,] WeylHeisenberg(double[] g0, int M)
        {
            int N = g0.Length, L = N / M;

            if (L <= 0)
                throw new Exception("Number of frequency shifts not defined correctly");

            Complex[,] G = new Complex[N, N];
            Complex c = 2 * Maths.Pi * Maths.I;
            double a = M / 2.0;

            Parallel.For(0, N, n =>
            {
                double phase = n - a / 2.0;
                int k, l, u, i, j;
                Complex exp, psi;

                for (k = 0; k < M; k++)
                {
                    exp = Maths.Exp(c * k / M * phase);

                    for (l = 0; l < L; l++)
                    {
                        u = l * M + k;
                        i = Maths.Mod(n - l * M, N);
                        j = Maths.Mod(n + M / 2 - l * M, N);

                        psi = new Complex(
                            (g0[i] * exp).Real,
                            (Maths.I * g0[j] * exp).Real);

                        G[n, u] = psi;
                    }
                }
            });

            return G;
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[] B = Matrice.Dot(A, U.Hermitian());
            return B;
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[] A = Matrice.Dot(B, U);
            return A;
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[,] V = WeylHeisenbergTransform.WeylHeisenberg(this.window, M, this.m, true);
            Complex[,] B;

            if (direction == Direction.Both)
            {
                B = U.Hermitian().Dot(A).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Hermitian().Dot(A);
            }
            else
            {
                B = A.Dot(V);
            }
            return B;
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[,] V = WeylHeisenbergTransform.WeylHeisenberg(this.window, M, this.m, true);
            Complex[,] A;

            if (direction == Direction.Both)
            {
                A = U.Dot(B).Dot(V.Hermitian());
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Dot(B);
            }
            else
            {
                A = B.Dot(V.Hermitian());
            }
            return A;
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Zak components
        private static FourierTransform DFT = new FourierTransform(false, Direction.Vertical);
        private static FastFourierTransform FFT = new FastFourierTransform(false, Direction.Vertical);

        /// <summary>
        /// Implements Zak-orthogonalization of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
        public static double[] Zak(double[] v, int M)
        {
            // Fast shaping orthogonalization algorithm
            // WH functions using a discrete Zak transform.
            // V.P. Volchkov, D.A. Petrov and V.M. Asiryan.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = v.Length;
            double[] vort = new double[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex[,] G = new Complex[L2, N];
            Complex[,] Z;

            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = v[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            if (Maths.IsPower(L2, 2))
            {
                Z = FFT.Forward(G);
            }
            else
            {
                Z = DFT.Forward(G);
            }

            double w = 2 / Math.Sqrt(M);
            double even, odd, phi;
            Complex z1, z2;

            for (i = 0; i < L; i++)
            {
                for (j = 0; j < N; j++)
                {
                    z1 = Z[i, j];
                    z2 = Z[L + i, j];

                    even = Math.Pow(z1.Abs, 2);
                    odd = Math.Pow(z2.Abs, 2);
                    phi = w / Math.Sqrt(even + odd);

                    Z[i, j] = z1 * phi;
                    Z[L + i, j] = z2 * phi;
                }
            }

            Complex sum;
            for (i = 0; i < N; i++)
            {
                sum = 0;
                for (j = 0; j < L2; j++)
                {
                    sum += Z[j, i];
                }
                vort[i] = (sum / L2).Real;
            }

            return vort;
        }
        /// <summary>
        /// Implements Zak-orthogonalization of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
        public static Complex[] Zak(Complex[] v, int M)
        {
            // Fast shaping orthogonalization algorithm
            // WH functions using a discrete Zak transform.
            // V.P. Volchkov, D.A. Petrov and V.M. Asiryan.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = v.Length;
            Complex[] vort = new Complex[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex[,] G = new Complex[L2, N];
            Complex[,] Z;

            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = v[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            if (Maths.IsPower(L2, 2))
            {
                Z = FFT.Forward(G);
            }
            else
            {
                Z = DFT.Forward(G);
            }

            double w = 2 / Math.Sqrt(M);
            double even, odd, phi;
            Complex z1, z2;

            for (i = 0; i < L; i++)
            {
                for (j = 0; j < N; j++)
                {
                    z1 = Z[i, j];
                    z2 = Z[L + i, j];

                    even = Math.Pow(z1.Abs, 2);
                    odd = Math.Pow(z2.Abs, 2);
                    phi = w / Math.Sqrt(even + odd);

                    Z[i, j] = z1 * phi;
                    Z[L + i, j] = z2 * phi;
                }
            }

            Complex sum;
            for (i = 0; i < N; i++)
            {
                sum = 0;
                for (j = 0; j < L2; j++)
                {
                    sum += Z[j, i];
                }
                vort[i] = sum / L2;
            }

            return vort;
        }
        #endregion
    }
}
