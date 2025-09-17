using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Discrete Chebyshev transform of the first kind (Type-I), orthonormal.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Chebyshev_transform
    /// </remarks>
    [Serializable]
    public class ChebyshevTransform : TransformBase, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the Chebyshev transform (Type-I, orthonormal).
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public ChebyshevTransform(Direction direction = Direction.Vertical)
        {
            this.Direction = direction;
        }
        #endregion

        #region Chebyshev static components
        /// <summary>
        /// Constructs the orthonormal Chebyshev (DCT-I) matrix of size n×n.
        /// </summary>
        /// <param name="n">Size (n ≥ 1)</param>
        /// <returns>Matrix H</returns>
        public static float[,] Matrix(int n)
        {
            if (n == 1)
                return new float[1, 1] { { 1f } };

            int M = n - 1;
            float[,] H = new float[n, n];
            float cm = Maths.Pi / M;                // π / M
            float norm = Maths.Sqrt(2.0f / M);      // √(2/M)
            float invSqrt2 = 1.0f / Maths.Sqrt2;    // 1/√2

            for (int k = 0; k < n; k++)
            {
                float ak = (k == 0 || k == M) ? invSqrt2 : 1.0f;
                for (int j = 0; j < n; j++)
                {
                    float aj = (j == 0 || j == M) ? invSqrt2 : 1.0f;
                    H[k, j] = norm * ak * aj * Maths.Cos(cm * k * j);
                }
            }
            return H;
        }
        #endregion

        #region Chebyshev Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            int N = A.Length;
            float[,] U = ChebyshevTransform.Matrix(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length;
            float[,] U = ChebyshevTransform.Matrix(N);
            return Matrice.Dot(B, Matrice.Transpose(U));
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = ChebyshevTransform.Matrix(N);
            float[,] V = ChebyshevTransform.Matrix(M);

            if (Direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transpose());
            }
            else if (Direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transpose());
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = ChebyshevTransform.Matrix(N);
            float[,] V = ChebyshevTransform.Matrix(M);

            if (Direction == Direction.Both)
            {
                return U.Transpose().Dot(B).Dot(V);
            }
            else if (Direction == Direction.Vertical)
            {
                return U.Transpose().Dot(B);
            }
            return B.Dot(V);
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            float[,] U = ChebyshevTransform.Matrix(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            float[,] U = ChebyshevTransform.Matrix(N);
            return Matrice.Dot(B, Matrice.Transpose(U));
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = ChebyshevTransform.Matrix(N);
            float[,] V = ChebyshevTransform.Matrix(M);

            if (Direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transpose());
            }
            else if (Direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transpose());
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = ChebyshevTransform.Matrix(N);
            float[,] V = ChebyshevTransform.Matrix(M);

            if (Direction == Direction.Both)
            {
                return U.Transpose().Dot(B).Dot(V);
            }
            else if (Direction == Direction.Vertical)
            {
                return U.Transpose().Dot(B);
            }
            return B.Dot(V);
        }
        #endregion
    }
}
