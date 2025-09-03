﻿using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Walsh-Hadamard transform.
    /// <remarks>
    /// More information can be found on the website:
    /// http://kibia.ru/teachers/kreindelin/pdf/2.pdf
    /// </remarks>
    /// </summary>
    [Serializable]
    public class WalshHadamardTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        private bool normalized;
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Walsh-Hadamard transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public WalshHadamardTransform(bool normalized = true, Direction direction = Direction.Vertical)
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

        #region Walsh-Hadamard static components
        /// <summary>
        /// Implements the construction of the Walsh-Hadamard matrix.
        /// </summary>
        /// <param name="powOf2">Power of 2</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int powOf2)
        {
            int iterations = powOf2 - 1;
            float[,] hadamard = WalshHadamardTransform.Matrix();

            if (iterations > 0)
            {
                float[,] H = WalshHadamardTransform.Matrix();

                for (int i = 0; i < iterations; i++)
                {
                    H = MatrixF.Kronecker(H, hadamard);
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
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            int N = A.Length;

            if (!MathF.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int n = (int)MathF.Log2(N);
            float[,] U = WalshHadamardTransform.Matrix(n);
            float[] B = MatrixF.Dot(A, U);

            if (normalized)
            {
                B = MatrixF.Div(B, MathF.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length;

            if (!MathF.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int n = (int)MathF.Log2(N);
            float[,] U = WalshHadamardTransform.Matrix(n);
            float[] A = MatrixF.Dot(B, U.Transponate());

            if (normalized)
            {
                A = MatrixF.Div(A, MathF.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);

            if (!MathF.IsPower(N, 2) || !MathF.IsPower(M, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[,] U = WalshHadamardTransform.Matrix((int)MathF.Log2(N));
            float[,] V = WalshHadamardTransform.Matrix((int)MathF.Log2(M));
            float[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transponate());
                B = normalized ? B.Div(MathF.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(MathF.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transponate());
                B = normalized ? B.Div(MathF.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);

            if (!MathF.IsPower(N, 2) || !MathF.IsPower(M, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[,] U = WalshHadamardTransform.Matrix((int)MathF.Log2(N));
            float[,] V = WalshHadamardTransform.Matrix((int)MathF.Log2(M));
            float[,] A;

            if (direction == Direction.Both)
            {
                A = U.Transponate().Dot(B).Dot(V);
                A = normalized ? A.Div(MathF.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Transponate().Dot(B);
                A = normalized ? A.Div(MathF.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = normalized ? A.Div(MathF.Sqrt(M)) : A;
            }

            return A;
        }
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public ComplexF[] Forward(ComplexF[] A)
        {
            int N = A.Length;

            if (!MathF.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int n = (int)MathF.Log2(N);
            float[,] U = WalshHadamardTransform.Matrix(n);
            ComplexF[] B = MatrixF.Dot(A, U);

            if (normalized)
            {
                B = MatrixF.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public ComplexF[] Backward(ComplexF[] B)
        {
            int N = B.Length;

            if (!MathF.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int n = (int)MathF.Log2(N);
            float[,] U = WalshHadamardTransform.Matrix(n);
            ComplexF[] A = MatrixF.Dot(B, U.Transponate());

            if (normalized)
            {
                A = MatrixF.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public ComplexF[,] Forward(ComplexF[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);

            if (!MathF.IsPower(N, 2) || !MathF.IsPower(M, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[,] U = WalshHadamardTransform.Matrix((int)MathF.Log2(N));
            float[,] V = WalshHadamardTransform.Matrix((int)MathF.Log2(M));
            ComplexF[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Math.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public ComplexF[,] Backward(ComplexF[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);

            if (!MathF.IsPower(N, 2) || !MathF.IsPower(M, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[,] U = WalshHadamardTransform.Matrix((int)MathF.Log2(N));
            float[,] V = WalshHadamardTransform.Matrix((int)MathF.Log2(M));
            ComplexF[,] A;

            if (direction == Direction.Both)
            {
                A = U.Transponate().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Transponate().Dot(B);
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
