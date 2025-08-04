using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Chebyshev transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Chebyshev_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ChebyshevTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Chebyshev transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public ChebyshevTransform(Direction direction = Direction.Vertical)
        {
            this.direction = direction;
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

        #region Static components
        /// <summary>
        /// Implements the construction of the Chebyshev transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int n)
        {
            float[,] T = new float[n, n];
            float[] x = ChebyshevNodes(n);

            for (int k = 0; k < n; k++)
            {
                float norm = k == 0 ? (float)Math.Sqrt(1.0 / n) : (float)Math.Sqrt(2.0 / n);

                for (int i = 0; i < n; i++)
                {
                    T[k, i] = norm * (float)Math.Cos(k * Math.Acos(x[i]));
                }
            }

            return T;
        }
        #endregion

        #region Chebyshev Transform
        /// <summary>
        /// Forward Chebyshev transform.
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
        /// Backward Chebyshev transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length;
            float[,] U = ChebyshevTransform.Matrix(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Forward Chebyshev transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = ChebyshevTransform.Matrix(N);
            float[,] V = ChebyshevTransform.Matrix(M);

            if (direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transponate());
        }
        /// <summary>
        /// Backward Chebyshev transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = ChebyshevTransform.Matrix(N);
            float[,] V = ChebyshevTransform.Matrix(M);

            if (direction == Direction.Both)
            {
                return U.Transponate().Dot(B).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                return U.Transponate().Dot(B);
            }
            return B.Dot(V);
        }
        /// <summary>
        /// Forward Chebyshev transform.
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
        /// Backward Chebyshev transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            float[,] U = ChebyshevTransform.Matrix(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Forward Chebyshev transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = ChebyshevTransform.Matrix(N);
            float[,] V = ChebyshevTransform.Matrix(M);

            if (direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transponate());
        }
        /// <summary>
        /// Backward Chebyshev transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = ChebyshevTransform.Matrix(N);
            float[,] V = ChebyshevTransform.Matrix(M);

            if (direction == Direction.Both)
            {
                return U.Transponate().Dot(B).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                return U.Transponate().Dot(B);
            }
            return B.Dot(V);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Returns the Chebyshev nodes of the first kind (orthogonal).
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        private static float[] ChebyshevNodes(int n)
        {
            float[] x = new float[n];

            for (int i = 0; i < n; i++)
            {
                x[i] = (float)Math.Cos(Math.PI * (i + 0.5) / n); // nodes on [-1, 1]
            }

            return x;
        }
        #endregion
    }
}
