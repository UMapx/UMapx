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
    /// https://ieeexplore.ieee.org/document/9117707/
    /// </remarks>
    /// </summary>
    [Serializable]
    public class WeylHeisenbergTransform : IWindowTransform, ITransform
    {
        #region Private data
        /// <summary>
        /// Windows function.
        /// </summary>
        protected IWindow window;
        /// <summary>
        /// Number of frequency shifts.
        /// </summary>
        protected int m;
        /// <summary>
        /// Processing direction.
        /// </summary>
        protected Direction direction;
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
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Matrix(float[] g0, int M)
        {
            int N = g0.Length, L = N / M;

            if (L <= 0)
                throw new Exception("Number of frequency shifts not defined correctly");

            Complex32[,] G = new Complex32[N, N * 2];
            Complex32 c = 2 * Maths.Pi * Maths.I;
            float a = M / 2.0f;

            Parallel.For(0, N, n =>
            {
                float phase = n - a / 2.0f;
                int k, l, u, i, j;
                Complex32 exp;

                for (k = 0; k < M; k++)
                {
                    exp = Maths.Exp(c * k / M * phase);

                    for (l = 0; l < L; l++)
                    {
                        u = l * M + k;
                        i = Maths.Mod(n - l * M, N);
                        j = Maths.Mod(n + M / 2 - l * M, N);

                        G[n, u + 0] =           g0[i] * exp / Maths.Sqrt(2);
                        G[n, u + N] = Maths.I * g0[j] * exp / Maths.Sqrt(2);
                    }
                }
            });

            return G;
        }
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
        public static Complex32[,] Matrix(IWindow window, int N, int M, bool orthogonalize = true)
        {
            return WeylHeisenbergTransform.Matrix(WeylHeisenbergTransform.GetPacket(window, N), M, orthogonalize);
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
        public static Complex32[,] Matrix(float[] g0, int M, bool orthogonalize = true)
        {
            if (orthogonalize)
            {
                ZakTransform zakTransform = new ZakTransform(M);

                return WeylHeisenbergTransform.Matrix(zakTransform.Orthogonalize(g0), M);
            }

            return WeylHeisenbergTransform.Matrix(g0, M);
        }
        /// <summary>
        /// Returns a vector of window function values.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="length">Number of samples</param>
        /// <returns>Array</returns>
        public static float[] GetPacket(IWindow window, int length)
        {
            // exeption by length
            if (window.FrameSize > length)
                return WeylHeisenbergTransform.NSymmetry(window, length);

            // params for approximation
            float[] w = WeylHeisenbergTransform.NSymmetry(window, (int)window.FrameSize);
            int n = w.Length;
            float min = Math.Min(w[0], w[n - 1]);
            float[] g = new float[length];
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
        private static float[] NSymmetry(IWindow window, int length)
        {
            // creaing window function
            float[] g = window.GetWindow(length + 1);
            float[] w = new float[length];

            // N-1 symmetric
            for (int i = 0; i < length; i++)
                w[i] = g[i];

            return w;
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public virtual Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            Complex32[,] U = WeylHeisenbergTransform.Matrix(this.window, N, this.m, true);
            Complex32[] B = Matrice.Dot(A, U.Hermitian());
            return B;
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public virtual Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            Complex32[,] U = WeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            Complex32[] A = Matrice.Dot(B, U);
            return A;
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public virtual Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            Complex32[,] U = WeylHeisenbergTransform.Matrix(this.window, N, this.m, true);
            Complex32[,] V = WeylHeisenbergTransform.Matrix(this.window, M, this.m, true);
            Complex32[,] B;

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
        public virtual Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            Complex32[,] U = WeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            Complex32[,] V = WeylHeisenbergTransform.Matrix(this.window, M / 2, this.m, true);
            Complex32[,] A;

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
        public virtual float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public virtual float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public virtual float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public virtual float[,] Backward(float[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
}
