using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines a group of real orthogonal bases and discrete Weyl-Heisenberg transforms.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://ieeexplore.ieee.org/document/8711969/
    /// </remarks>
    [Serializable]
    public class RealWeylHeisenbergTransform : WeylHeisenbergTransform, IWindowTransform, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes a group of real orthogonal bases and Weyl-Heisenberg transformations.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N/4]</param>
        /// <param name="direction">Processing direction</param>
        public RealWeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical) : base(window, m, direction) { }
        #endregion

        #region Weyl-Heisenberg static components
        /// <summary>
        /// Returns the real Weyl-Heisenberg basis matrix.
        /// </summary>
        /// <remarks>
        /// Matrix dimension [2N, 2N], where N = M * L.
        /// </remarks>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts [4, N/4]</param>
        /// <returns>Matrix</returns>
        public new static float[,] Matrix(float[] g0, int M)
        {
            int N = g0.Length, L = N / M;

            if (L <= 0)
                throw new ArgumentException("Number of frequency shifts not defined correctly");

            float[,] G = new float[2 * N, 2 * N];
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

                        var G1 =           g0[i] * exp;
                        var G2 = Maths.I * g0[j] * exp;

                        // implements the [2N, 2N] WH matrix
                        // https://github.com/asiryan/Weyl-Heisenberg-Toolbox/blob/master/matlab/toolbox_scripts/weylhzr.m

                        G[n, u + 0] = G1.Real; G[n + N, u + 0] = G1.Imag;
                        G[n, u + N] = G2.Real; G[n + N, u + N] = G2.Imag;
                    }
                }
            });

            return G;
        }
        /// <summary>
        /// Returns the real Weyl-Heisenberg basis matrix.
        /// </summary>
        /// <remarks>
        /// Matrix dimension [2N, 2N], where N = M * L.
        /// </remarks>
        /// <param name="window">Windows function</param>
        /// <param name="N">Number of samples</param>
        /// <param name="M">Number of frequency shifts [4, N/4]</param>
        /// <param name="orthogonalize">Orthogonalized matrix or not</param>
        /// <returns>Matrix</returns>
        public new static float[,] Matrix(IWindow window, int N, int M, bool orthogonalize = true)
        {
            return RealWeylHeisenbergTransform.Matrix(WeylHeisenbergTransform.Packet(window, N), M, orthogonalize);
        }
        /// <summary>
        /// Returns the real Weyl-Heisenberg basis matrix.
        /// </summary>
        /// <remarks>
        /// Matrix dimension [2N, 2N], where N = M * L.
        /// </remarks>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts [4, N/4]</param>
        /// <param name="orthogonalize">Orthogonalized matrix or not</param>
        /// <returns>Matrix</returns>
        public new static float[,] Matrix(float[] g0, int M, bool orthogonalize = true)
        {
            if (orthogonalize)
            {
                ZakTransform zakTransform = new ZakTransform(M);

                return RealWeylHeisenbergTransform.Matrix(zakTransform.Orthogonalize(g0), M);
            }

            return RealWeylHeisenbergTransform.Matrix(g0, M);
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override float[] Forward(float[] A)
        {
            int N = A.Length;
            float[,] U = RealWeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            float[] B = Core.Matrice.Dot(A, (float[,])Core.Matrice.Transponate(U));
            return B;
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override float[] Backward(float[] B)
        {
            int N = B.Length;
            float[,] U = RealWeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            float[] A = Core.Matrice.Dot(B, U);
            return A;
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public override float[,] Forward(float[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = RealWeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            float[,] V = RealWeylHeisenbergTransform.Matrix(this.window, M / 2, this.m, true);
            float[,] B;

            if (direction == Direction.Both)
            {
                B = U.Transponate().Dot(A).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Transponate().Dot(A);
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
        public override float[,] Backward(float[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = RealWeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            float[,] V = RealWeylHeisenbergTransform.Matrix(this.window, M / 2, this.m, true);
            float[,] A;

            if (direction == Direction.Both)
            {
                A = U.Dot(B).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Dot(B);
            }
            else
            {
                A = B.Dot(V.Transponate());
            }
            return A.Div(2);
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            float[,] U = RealWeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            Complex32[] B = Core.Matrice.Dot(A, (float[,])Core.Matrice.Transponate(U));
            return B;
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            float[,] U = RealWeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            Complex32[] A = Core.Matrice.Dot(B, U);
            return A;
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public override Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = RealWeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            float[,] V = RealWeylHeisenbergTransform.Matrix(this.window, M / 2, this.m, true);
            Complex32[,] B;

            if (direction == Direction.Both)
            {
                B = U.Transponate().Dot(A).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Transponate().Dot(A);
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
        public override Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = RealWeylHeisenbergTransform.Matrix(this.window, N / 2, this.m, true);
            float[,] V = RealWeylHeisenbergTransform.Matrix(this.window, M / 2, this.m, true);
            Complex32[,] A;

            if (direction == Direction.Both)
            {
                A = U.Dot(B).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Dot(B);
            }
            else
            {
                A = B.Dot(V.Transponate());
            }
            return A;
        }
        #endregion
    }
}
