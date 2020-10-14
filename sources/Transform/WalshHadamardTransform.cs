using System;
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
        public static double[,] Hadamard(int powOf2)
        {
            int iterations = powOf2 - 1;
            double[,] hadamard = WalshHadamardTransform.Hadamard();

            if (iterations > 0)
            {
                double[,] H = WalshHadamardTransform.Hadamard();

                for (int i = 0; i < iterations; i++)
                {
                    H = Matrice.Kronecker(H, hadamard);
                }
                return H;
            }

            return hadamard;
        }
        /// <summary>
        /// Implements the construction of the Walsh-Hadamard matrix [2 x 2].
        /// </summary>
        /// <returns>Matrix</returns>
        public static double[,] Hadamard()
        {
            return (new double[2, 2] { { 1, 1 }, { 1, -1 } });
        }
        #endregion

        #region Walsh-Hadamard Transform
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length;

            if (!Maths.IsPower(N, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            int n = (int)Maths.Log2(N);
            double[,] U = WalshHadamardTransform.Hadamard(n);
            double[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length;

            if (!Maths.IsPower(N, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            int n = (int)Maths.Log2(N);
            double[,] U = WalshHadamardTransform.Hadamard(n);
            double[] A = Matrice.Dot(B, U.Transponate());

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            double[,] U = WalshHadamardTransform.Hadamard((int)Maths.Log2(N));
            double[,] V = WalshHadamardTransform.Hadamard((int)Maths.Log2(M));
            double[,] B;

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
        public double[,] Backward(double[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            double[,] U = WalshHadamardTransform.Hadamard((int)Maths.Log2(N));
            double[,] V = WalshHadamardTransform.Hadamard((int)Maths.Log2(M));
            double[,] A;

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
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;

            if (!Maths.IsPower(N, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            int n = (int)Maths.Log2(N);
            double[,] U = WalshHadamardTransform.Hadamard(n);
            Complex[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Backward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;

            if (!Maths.IsPower(N, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            int n = (int)Maths.Log2(N);
            double[,] U = WalshHadamardTransform.Hadamard(n);
            Complex[] A = Matrice.Dot(B, U.Transponate());

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Forward Walsh-Hadamard transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            double[,] U = WalshHadamardTransform.Hadamard((int)Maths.Log2(N));
            double[,] V = WalshHadamardTransform.Hadamard((int)Maths.Log2(M));
            Complex[,] B;

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
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new Exception("Dimension of the signal must be a power of 2");

            double[,] U = WalshHadamardTransform.Hadamard((int)Maths.Log2(N));
            double[,] V = WalshHadamardTransform.Hadamard((int)Maths.Log2(M));
            Complex[,] A;

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
