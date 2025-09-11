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
    public class RealWeylHeisenbergHartleyTransform : WeylHeisenbergTransform, IWindowTransform, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes a group of real orthogonal bases and Weyl-Heisenberg transformations.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N/4]</param>
        /// <param name="direction">Processing direction</param>
        public RealWeylHeisenbergHartleyTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical) : base(window, m, direction) { }
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
            int N = g0.Length;
            if (M <= 0 || N % M != 0)
                throw new ArgumentException("Number of frequency shifts not defined correctly");
            int L = N / M;

            // Optional: for exact phasing used below, M should be even; for WH real blocks typically M % 4 == 0.
            if ((M & 1) != 0)
                throw new ArgumentException("M must be even (preferably divisible by 4) for this phasing.");

            float[,] G = new float[2 * N, 2 * N];

            // Constants
            float twoPiOverM = 2.0f * Maths.Pi / M;
            float halfPi = 0.5f * Maths.Pi;
            float quarterM = 0.25f * M; // shift of M/4 used in 'phase'

            Parallel.For(0, N, n =>
            {
                // phase = n - M/4  (same as original: a = M/2, then phase = n - a/2)
                float phase = n - quarterM;

                for (int k = 0; k < M; k++)
                {
                    // θ = 2π * (k/M) * phase
                    float theta = twoPiOverM * k * phase;

                    // Hartley evaluations
                    float cas1 = Special.Cas(theta);             // = cosθ + sinθ
                    float cas2 = Special.Cas(theta - halfPi);    // = sinθ - cosθ

                    // Recover cosθ and sinθ from cas:
                    float c = 0.5f * (cas1 - cas2); // cosθ
                    float s = 0.5f * (cas1 + cas2); // sinθ

                    for (int l = 0; l < L; l++)
                    {
                        int u = l * M + k;
                        int i = Maths.Mod(n - l * M, N);
                        int j = Maths.Mod(n + (M / 2) - l * M, N);

                        float gi = g0[i];
                        float gj = g0[j];

                        // This exactly matches:
                        // G1 = gi * (cosθ + i sinθ)
                        // G2 = i * gj * (cosθ + i sinθ) = gj * (cos(θ+π/2) + i sin(θ+π/2))
                        // => Re(G2) = -gj * sinθ, Im(G2) = gj * cosθ
                        G[n, u + 0] = gi * c;   // Re(G1)
                        G[n + N, u + 0] = gi * s;   // Im(G1)

                        G[n, u + N] = -gj * s;  // Re(G2)
                        G[n + N, u + N] = gj * c;  // Im(G2)
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
            float[] B = Core.Matrice.Dot(A, (float[,])Core.Matrice.Transpose(U));
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

            if (Direction == Direction.Both)
            {
                B = U.Transpose().Dot(A).Dot(V);
            }
            else if (Direction == Direction.Vertical)
            {
                B = U.Transpose().Dot(A);
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

            if (Direction == Direction.Both)
            {
                A = U.Dot(B).Dot(V.Transpose());
            }
            else if (Direction == Direction.Vertical)
            {
                A = U.Dot(B);
            }
            else
            {
                A = B.Dot(V.Transpose());
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
            Complex32[] B = Core.Matrice.Dot(A, (float[,])Core.Matrice.Transpose(U));
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

            if (Direction == Direction.Both)
            {
                B = U.Transpose().Dot(A).Dot(V);
            }
            else if (Direction == Direction.Vertical)
            {
                B = U.Transpose().Dot(A);
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

            if (Direction == Direction.Both)
            {
                A = U.Dot(B).Dot(V.Transpose());
            }
            else if (Direction == Direction.Vertical)
            {
                A = U.Dot(B);
            }
            else
            {
                A = B.Dot(V.Transpose());
            }
            return A;
        }
        #endregion
    }
}
