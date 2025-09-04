using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Fourier transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Fourier_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FourierTransform : ITransform
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
        /// Initializes the Fourier transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FourierTransform(bool normalized = true, Direction direction = Direction.Vertical)
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

        #region Fourier static components
        /// <summary>
        /// Implements the construction of the Fourier matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Matrix(int n)
        {
            Complex32[,] H = new Complex32[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Maths.Exp(-2 * Maths.Pi * Maths.I * i * j / n);
                }
            }
            return H;
        }
        #endregion

        #region Fourier Transform
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            Complex32[,] U = FourierTransform.Matrix(N);
            Complex32[] B = Core.Matrix.Dot(A, U);

            if (normalized)
            {
                B = Core.Matrix.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            Complex32[,] U = FourierTransform.Matrix(N);
            Complex32[] A = Core.Matrix.Dot(B, Core.Matrix.Hermitian(U));

            if (normalized)
            {
                A = Core.Matrix.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            Complex32[,] U = FourierTransform.Matrix(N);
            Complex32[,] V = FourierTransform.Matrix(M);
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
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            Complex32[,] U = FourierTransform.Matrix(N);
            Complex32[,] V = FourierTransform.Matrix(M);
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
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Fourier transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Fourier transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Fourier transform.
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
