using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the sine transform.
    /// <remarks>
    /// More information can be found on the website:
    /// http://sernam.ru/book_prett1.php?id=91
    /// </remarks>
    /// </summary>
    [Serializable]
    public class SineTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the sine transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public SineTransform(Direction direction = Direction.Vertical)
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

        #region Sine static components
        /// <summary>
        /// Implements the construction of the sine transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static double[,] Matrix(int n)
        {
            int j, i;
            double[,] H = new double[n, n];
            double n1 = n + 1;
            double scale = Maths.Sqrt(2.0 / n1);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Math.Sin(Maths.Pi * (j + 1) * (i + 1) / n1) * scale;
                }
            }

            return (H);
        }
        #endregion

        #region Sine Transform
        /// <summary>
        /// Forward sine transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length;
            double[,] U = SineTransform.Matrix(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward sine transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length;
            double[,] U = SineTransform.Matrix(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Forward sine transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = SineTransform.Matrix(N);
            double[,] V = SineTransform.Matrix(M);

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
        /// Backward sine transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = SineTransform.Matrix(N);
            double[,] V = SineTransform.Matrix(M);

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
        /// Forward sine transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            double[,] U = SineTransform.Matrix(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward sine transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            double[,] U = SineTransform.Matrix(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Forward sine transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = SineTransform.Matrix(N);
            double[,] V = SineTransform.Matrix(M);

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
        /// Backward sine transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = SineTransform.Matrix(N);
            double[,] V = SineTransform.Matrix(M);

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
