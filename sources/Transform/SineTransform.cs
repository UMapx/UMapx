﻿using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the sine transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://sernam.ru/book_prett1.php?id=91
    /// </remarks>
    [Serializable]
    public class SineTransform : TransformBase, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the sine transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public SineTransform(Direction direction = Direction.Vertical)
        {
            this.Direction = direction;
        }
        #endregion

        #region Sine static components
        /// <summary>
        /// Implements the construction of the sine transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int n)
        {
            int j, i;
            float[,] H = new float[n, n];
            float n1 = n + 1;
            float scale = Maths.Sqrt(2.0f / n1);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Maths.Sin(Maths.Pi * (j + 1) * (i + 1) / n1) * scale;
                }
            }

            return H;
        }
        #endregion

        #region Sine Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            int N = A.Length;
            float[,] U = SineTransform.Matrix(N);
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
            float[,] U = SineTransform.Matrix(N);
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
            float[,] U = SineTransform.Matrix(N);
            float[,] V = SineTransform.Matrix(M);

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
            float[,] U = SineTransform.Matrix(N);
            float[,] V = SineTransform.Matrix(M);

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
            float[,] U = SineTransform.Matrix(N);
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
            float[,] U = SineTransform.Matrix(N);
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
            float[,] U = SineTransform.Matrix(N);
            float[,] V = SineTransform.Matrix(M);

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
            float[,] U = SineTransform.Matrix(N);
            float[,] V = SineTransform.Matrix(M);

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
