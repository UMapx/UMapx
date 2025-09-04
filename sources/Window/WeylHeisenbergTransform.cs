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
    /// https://ieeexplore.ieee.org/document/8711969/
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
        /// <param name="m">Number of frequency shifts [4, N/2]</param>
        /// <param name="direction">Processing direction</param>
        public WeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical)
        {
            Window = window; M = m; Direction = direction;
        }
        /// <summary>
        /// Gets or sets number of frequency shifts [4, N].
        /// </summary>
        public int M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value <= 2)
                    throw new ArgumentException("M must be greater than 2");

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
        /// Matrix dimension [N, 2N], where N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Matrix</returns>
        public static ComplexF[,] Matrix(float[] g0, int M)
        {
            int N = g0.Length, L = N / M;

            if (L <= 0)
                throw new ArgumentException("Number of frequency shifts not defined correctly");

            ComplexF[,] G = new ComplexF[N, 2 * N];
            ComplexF c = 2 * MathsF.Pi * MathsF.I;
            float a = M / 2.0f;

            Parallel.For(0, N, n =>
            {
                float phase = n - a / 2.0f;
                int k, l, u, i, j;
                ComplexF exp;

                for (k = 0; k < M; k++)
                {
                    exp = MathsF.Exp(c * k / M * phase);

                    for (l = 0; l < L; l++)
                    {
                        u = l * M + k;
                        i = MathsF.Mod(n - l * M, N);
                        j = MathsF.Mod(n + M / 2 - l * M, N);

                        // implements the [N, 2N] WH matrix
                        // https://github.com/asiryan/Weyl-Heisenberg-Toolbox/blob/master/matlab/toolbox_scripts/weylhzf.m

                        G[n, u + 0] =           g0[i] * exp;
                        G[n, u + N] = MathsF.I * g0[j] * exp;
                    }
                }
            });

            return G;
        }
        /// <summary>
        /// Returns the complex Weyl-Heisenberg basis matrix.
        /// <remarks>
        /// Matrix dimension [N, 2N], where N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="N">Number of samples</param>
        /// <param name="M">Number of frequency shifts [4, N/2]</param>
        /// <param name="orthogonalize">Orthogonalized matrix or not</param>
        /// <returns>Matrix</returns>
        public static ComplexF[,] Matrix(IWindow window, int N, int M, bool orthogonalize = true)
        {
            return WeylHeisenbergTransform.Matrix(WeylHeisenbergTransform.Packet(window, N), M, orthogonalize);
        }
        /// <summary>
        /// Returns the complex Weyl-Heisenberg basis matrix.
        /// <remarks>
        /// Matrix dimension [N, 2N], where N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts [4, N/2]</param>
        /// <param name="orthogonalize">Orthogonalized matrix or not</param>
        /// <returns>Matrix</returns>
        public static ComplexF[,] Matrix(float[] g0, int M, bool orthogonalize = true)
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
        public static float[] Packet(IWindow window, int length)
        {
            // exception by length
            if (window.FrameSize >= length)
                return window.GetWindow(length);

            var baseLen = window.FrameSize;
            var g = window.GetWindow(baseLen);

            // Target buffer and symmetric padding sizes
            int pad = length - baseLen;
            int leftPad = pad / 2;
            int rightPad = pad - leftPad;

            var w = new float[length];

            // Use the window's own edge levels as "its minima" for padding
            // (for most windows these are indeed the minima).
            float leftVal = g[0];
            float rightVal = g[baseLen - 1];

            // Left pad
            for (int i = 0; i < leftPad; i++)
                w[i] = leftVal;

            // Copy the core window centered
            Buffer.BlockCopy(g, 0, w, sizeof(float) * leftPad, baseLen * sizeof(float));

            // Right pad
            for (int i = length - rightPad; i < length; i++)
                w[i] = rightVal;

            return w;
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public virtual ComplexF[] Forward(ComplexF[] A)
        {
            int N = A.Length;
            ComplexF[,] U = WeylHeisenbergTransform.Matrix(this.window, N, this.m, true);
            ComplexF[] B = MatrixF.Dot(A, U.Hermitian());
            return B;
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public virtual ComplexF[] Backward(ComplexF[] B)
        {
            int N = B.Length;
            ComplexF[,] U = WeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            ComplexF[] A = MatrixF.Dot(B, U);
            return A.Div(2);
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public virtual ComplexF[,] Forward(ComplexF[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            ComplexF[,] U = WeylHeisenbergTransform.Matrix(this.window, N, this.m, true);
            ComplexF[,] V = WeylHeisenbergTransform.Matrix(this.window, M, this.m, true);
            ComplexF[,] B;

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
        public virtual ComplexF[,] Backward(ComplexF[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            ComplexF[,] U = WeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            ComplexF[,] V = WeylHeisenbergTransform.Matrix(this.window, M / 2, this.m, true);
            ComplexF[,] A;

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
            return A.Div(2);
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
