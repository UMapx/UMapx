using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Hartley transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastHartleyTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private FastFourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Hartley transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastHartleyTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(normalized, direction);
        }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        /// <summary>
        /// Gets or sets the processing direction.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.FFT.Direction;
            }
            set
            {
                this.FFT.Direction = value;
            }
        }
        #endregion

        #region Fast Hartley Transform
        /// <summary>
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            Complex[] B = Matrice.ToComplex(A);
            B = FFT.Forward(B);

            int length = A.Length, i;
            double[] Hk = new double[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = B[i].Real + B[i].Imag;
            }

            return Hk;
        }
        /// <summary>
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            Complex[] A = Matrice.ToComplex(B);
            A = FFT.Backward(A);

            int length = B.Length, i;
            double[] Hk = new double[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = A[i].Real - A[i].Imag;
            }

            return Hk;
        }
        /// <summary>
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            Complex[,] B = Matrice.ToComplex(A);
            B = FFT.Forward(B);

            int width = A.GetLength(1), height = A.GetLength(0);
            double[,] Hk = new double[height, width];
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    Hk[i, j] = B[i, j].Real + B[i, j].Imag;
                }
            }

            return Hk;
        }
        /// <summary>
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            Complex[,] A = Matrice.ToComplex(B);
            A = FFT.Backward(A);

            int width = B.GetLength(1), height = B.GetLength(0);
            double[,] Hk = new double[height, width];
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    Hk[i, j] = A[i, j].Real - A[i, j].Imag;
                }
            }

            return Hk;
        }
        /// <summary>
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            int length = A.Length, i;
            Complex[] B = FFT.Forward(A);
            Complex[] Hk = new Complex[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = B[i].Real + B[i].Imag;
            }

            return Hk;
        }
        /// <summary>
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            int length = B.Length, i;
            Complex[] A = FFT.Backward(B);
            Complex[] Hk = new Complex[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = A[i].Real - A[i].Imag;
            }

            return Hk;
        }
        /// <summary>
        /// Forward Hartley transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = FFT.Forward(A);
            int width = A.GetLength(1), height = A.GetLength(0);
            Complex[,] Hk = new Complex[height, width];
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    Hk[i, j] = B[i, j].Real + B[i, j].Imag;
                }
            }

            return Hk;
        }
        /// <summary>
        /// Backward Hartley transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = FFT.Backward(B);
            int width = B.GetLength(1), height = B.GetLength(0);
            Complex[,] Hk = new Complex[height, width];
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    Hk[i, j] = A[i, j].Real - A[i, j].Imag;
                }
            }

            return Hk;
        }
        #endregion
    }
}
