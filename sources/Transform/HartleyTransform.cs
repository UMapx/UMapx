﻿using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Hartley transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
    /// </remarks>
    [Serializable]
    public class HartleyTransform : ITransform
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
        /// Initializes the Hartley transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public HartleyTransform(bool normalized = true, Direction direction = Direction.Vertical)
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

        #region Hartley static components
        /// <summary>
        /// Implements the construction of the Hartley transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int n)
        {
            int j, i;
            float[,] H = new float[n, n];
            float s = (2.0f * Maths.Pi) / n;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Special.Cas(s * i * j);
                }
            }

            return (H);
        }
        #endregion

        #region Hartley Transform
        /// <summary>
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            int N = A.Length;
            float[,] U = HartleyTransform.Matrix(N);
            float[] B = Core.Matrice.Dot(A, U);

            if (normalized)
            {
                B = Core.Matrice.Div(B, Maths.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length;
            float[,] U = HartleyTransform.Matrix(N);
            float[] A = Core.Matrice.Dot(B, (float[,])Core.Matrice.Transpose(U));

            if (normalized)
            {
                A = Core.Matrice.Div(A, Maths.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = HartleyTransform.Matrix(N);
            float[,] V = HartleyTransform.Matrix(M);
            float[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transpose());
                B = normalized ? B.Div(Maths.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
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
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = HartleyTransform.Matrix(N);
            float[,] V = HartleyTransform.Matrix(M);
            float[,] A;

            if (direction == Direction.Both)
            {
                A = U.Transpose().Dot(B).Dot(V);
                A = normalized ? A.Div(Maths.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
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
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            float[,] U = HartleyTransform.Matrix(N);
            Complex32[] B = Core.Matrice.Dot(A, U);

            if (normalized)
            {
                B = Core.Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            float[,] U = HartleyTransform.Matrix(N);
            Complex32[] A = Core.Matrice.Dot(B, (float[,])Core.Matrice.Transpose(U));

            if (normalized)
            {
                A = Core.Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = HartleyTransform.Matrix(N);
            float[,] V = HartleyTransform.Matrix(M);
            Complex32[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transpose());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
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
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = HartleyTransform.Matrix(N);
            float[,] V = HartleyTransform.Matrix(M);
            Complex32[,] A;

            if (direction == Direction.Both)
            {
                A = U.Transpose().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
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
