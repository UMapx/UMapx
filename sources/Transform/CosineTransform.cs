using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the cosine transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_cosine_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class CosineTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the cosine transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public CosineTransform(Direction direction = Direction.Vertical)
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

        #region Cosine static components
        /// <summary>
        /// Implements the construction of the cosine transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static double[,] Matrix(int n)
        {
            int j, i;
            double[,] H = new double[n, n];
            double c = Maths.Pi / (2.0 * n);
            double g1 = Maths.Sqrt(1.0 / n);
            double g2 = Math.Sqrt(2.0 / n);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = (i == 0) ? g1 : Math.Cos((2 * j + 1) * i * c) * g2;
                }
            }

            return (H);
        }
        #endregion

        #region Cosine Transform
        /// <summary>
        /// Forward cosine transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length;
            double[,] U = CosineTransform.Matrix(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward cosine transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length;
            double[,] U = CosineTransform.Matrix(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Forward cosine transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = CosineTransform.Matrix(N);
            double[,] V = CosineTransform.Matrix(M);

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
        /// Backward cosine transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = CosineTransform.Matrix(N);
            double[,] V = CosineTransform.Matrix(M);

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
        /// Forward cosine transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            double[,] U = CosineTransform.Matrix(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward cosine transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            double[,] U = CosineTransform.Matrix(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Forward cosine transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = CosineTransform.Matrix(N);
            double[,] V = CosineTransform.Matrix(M);

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
        /// Backward cosine transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = CosineTransform.Matrix(N);
            double[,] V = CosineTransform.Matrix(M);

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
    }
}
