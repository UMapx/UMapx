﻿using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Laplace transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Laplace_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LaplaceTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Standard deviation.
        /// </summary>
        private float sigma;
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
        /// Initializes the Laplace transform.
        /// </summary>
        /// <param name="sigma">Standard deviation (0, 1)</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public LaplaceTransform(float sigma = 0.0005f, bool normalized = true, Direction direction = Direction.Vertical)
        {
            Sigma = sigma; this.normalized = normalized; this.direction = direction;
        }
        /// <summary>
        /// Gets or sets the standard deviation (0, 1).
        /// <remarks>
        /// If σ = 0, then the Laplace transform takes the form of a Fourier transform.
        /// </remarks>
        /// </summary>
        public float Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Invalid argument value");

                this.sigma = value;
            }
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

        #region Laplace static components
        /// <summary>
        /// Implements the construction of the Laplace matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <param name="sigma">Standard deviation (0, 1)</param>
        /// <param name="backward">Return backward transformation matrix or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Matrix(int n, float sigma, bool backward = false)
        {
            Complex32[,] H = new Complex32[n, n];
            float factor;
            int i, j;

            // inverse matrix or not?
            if (backward)
            {
                for (i = 0; i < n; i++)
                {
                    factor = Maths.Exp(-sigma * i);

                    for (j = 0; j < n; j++)
                    {
                        H[i, j] = Maths.Exp(-2 * Maths.Pi * Maths.I * j / n * i) / factor;
                    }
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    factor = Maths.Exp(-sigma * i);

                    for (j = 0; j < n; j++)
                    {
                        H[i, j] = factor * Maths.Exp(-2 * Maths.Pi * Maths.I * j / n * i);
                    }
                }
            }
            return H;
        }
        #endregion

        #region Laplace Transform
        /// <summary>
        /// Forward Laplace transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            Complex32[,] U = LaplaceTransform.Matrix(N, sigma);
            Complex32[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Laplace transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            Complex32[,] U = LaplaceTransform.Matrix(N, sigma, true);
            Complex32[] A = Matrice.Dot(B, U.Hermitian());

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Laplace transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            Complex32[,] U = LaplaceTransform.Matrix(N, sigma);
            Complex32[,] V = LaplaceTransform.Matrix(M, sigma);
            Complex32[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Hermitian());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Math.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Hermitian());
                B = normalized ? B.Div(Math.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Backward Laplace transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            Complex32[,] U = LaplaceTransform.Matrix(N, sigma, true);
            Complex32[,] V = LaplaceTransform.Matrix(M, sigma, true);
            Complex32[,] A;

            if (direction == Direction.Both)
            {
                A = U.Hermitian().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Hermitian().Dot(B);
                A = normalized ? A.Div(Math.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = normalized ? A.Div(Math.Sqrt(M)) : A;
            }

            return A;
        }
        /// <summary>
        /// Forward Laplace transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Laplace transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Laplace transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Laplace transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
}
