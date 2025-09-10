using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Walsh-Hadamard transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://kibia.ru/teachers/kreindelin/pdf/2.pdf
    /// </remarks>
    [Serializable]
    public class WalshHadamardTransform : TransformBase, ITransform
    {
        #region Private data
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        private bool normalized;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Walsh-Hadamard transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public WalshHadamardTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.normalized = normalized; this.Direction = direction;
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
        #endregion

        #region Walsh-Hadamard static components
        /// <summary>
        /// Implements the construction of the Walsh-Hadamard matrix.
        /// </summary>
        /// <param name="powOf2">Power of 2</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int powOf2)
        {
            if (powOf2 < 0)
                throw new ArgumentOutOfRangeException(nameof(powOf2), "Power must be non-negative");

            if (powOf2 == 0)
                return new float[1, 1] { { 1 } };

            int iterations = powOf2 - 1;
            float[,] hadamard = WalshHadamardTransform.Matrix();

            if (iterations > 0)
            {
                float[,] H = WalshHadamardTransform.Matrix();

                for (int i = 0; i < iterations; i++)
                {
                    H = Core.Matrice.Kronecker(H, hadamard);
                }
                return H;
            }

            return hadamard;
        }
        /// <summary>
        /// Implements the construction of the Walsh-Hadamard matrix [2 x 2].
        /// </summary>
        /// <returns>Matrix</returns>
        public static float[,] Matrix()
        {
            return new float[2, 2] { { 1, 1 }, { 1, -1 } };
        }
        #endregion

        #region Walsh-Hadamard Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            int N = A.Length;

            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int n = (int)Maths.Log2(N);
            float[,] U = WalshHadamardTransform.Matrix(n);
            float[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Maths.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length;

            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int n = (int)Maths.Log2(N);
            float[,] U = WalshHadamardTransform.Matrix(n);
            float[] A = Matrice.Dot(B, Matrice.Transpose(U));

            if (normalized)
            {
                A = Matrice.Div(A, Maths.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[,] U = WalshHadamardTransform.Matrix((int)Maths.Log2(N));
            float[,] V = WalshHadamardTransform.Matrix((int)Maths.Log2(M));
            float[,] B;

            if (Direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transpose());
                B = normalized ? B.Div(Maths.Sqrt(N * M)) : B;
            }
            else if (Direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Maths.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transpose());
                B = normalized ? B.Div(Maths.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[,] U = WalshHadamardTransform.Matrix((int)Maths.Log2(N));
            float[,] V = WalshHadamardTransform.Matrix((int)Maths.Log2(M));
            float[,] A;

            if (Direction == Direction.Both)
            {
                A = U.Transpose().Dot(B).Dot(V);
                A = normalized ? A.Div(Maths.Sqrt(N * M)) : A;
            }
            else if (Direction == Direction.Vertical)
            {
                A = U.Transpose().Dot(B);
                A = normalized ? A.Div(Maths.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = normalized ? A.Div(Maths.Sqrt(M)) : A;
            }

            return A;
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;

            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int n = (int)Maths.Log2(N);
            float[,] U = WalshHadamardTransform.Matrix(n);
            Complex32[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;

            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int n = (int)Maths.Log2(N);
            float[,] U = WalshHadamardTransform.Matrix(n);
            Complex32[] A = Matrice.Dot(B, Matrice.Transpose(U));

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[,] U = WalshHadamardTransform.Matrix((int)Maths.Log2(N));
            float[,] V = WalshHadamardTransform.Matrix((int)Maths.Log2(M));
            Complex32[,] B;

            if (Direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transpose());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (Direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Math.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transpose());
                B = normalized ? B.Div(Math.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[,] U = WalshHadamardTransform.Matrix((int)Maths.Log2(N));
            float[,] V = WalshHadamardTransform.Matrix((int)Maths.Log2(M));
            Complex32[,] A;

            if (Direction == Direction.Both)
            {
                A = U.Transpose().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (Direction == Direction.Vertical)
            {
                A = U.Transpose().Dot(B);
                A = normalized ? A.Div(Math.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = normalized ? A.Div(Math.Sqrt(M)) : A;
            }

            return A;
        }
        #endregion
    }
}
