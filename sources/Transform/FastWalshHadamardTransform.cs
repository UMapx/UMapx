using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Walsh-Hadamard transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.mathworks.com/matlabcentral/fileexchange/6879-fast-walsh-hadamard-transform
    /// </remarks>
    [Serializable]
    public class FastWalshHadamardTransform : FloatToComplex32TransformBase, ITransform
    {
        #region Private data
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        private bool normalized;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Walsh-Hadamard transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastWalshHadamardTransform(bool normalized = true, Direction direction = Direction.Vertical)
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

        #region Fast Walsh-Hadamard Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override float[] Forward(float[] A)
        {
            int N = A.Length;
            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[] B = (float[])A.Clone(); FWHT(B);

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
        public override float[] Backward(float[] B)
        {
            int N = B.Length;
            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            float[] A = (float[])B.Clone(); FWHT(A);

            if (normalized)
            {
                A = Matrice.Div(A, Maths.Sqrt(N));
            }

            return A;
        }
        #endregion

        #region Private data
        /// <summary>
        /// Fast Walsh-Hadamard transform.
        /// </summary>
        /// <param name="data">Array</param>
        private void FWHT(float[] data)
        {
            int N = data.Length;

            if (N == 1)
                return;

            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int log2N = (int)Maths.Log2(N);
            float x_even, x_odd;

            int k0 = N, k1 = 1, k2 = N / 2;
            int x, y, z, i, j, l;

            for (x = 0; x < log2N; x++)
            {
                l = 0;

                for (y = 0; y < k1; y++, l += k0)
                {
                    for (z = 0; z < k2; z++)
                    {
                        i = z + l; j = i + k2;

                        x_even = data[i];
                        x_odd = data[j];

                        data[i] = x_even + x_odd;
                        data[j] = x_even - x_odd;
                    }
                }

                k0 /= 2; k1 *= 2; k2 /= 2;
            }
        }
        /// <summary>
        /// Fast Walsh-Hadamard transform.
        /// </summary>
        /// <param name="data">Array</param>
        private void FWHT(Complex32[] data)
        {
            int N = data.Length;

            if (N == 1)
                return;

            if (!Maths.IsPower(N, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int log2N = (int)Maths.Log2(N);
            Complex32 x_even, x_odd;

            int k0 = N, k1 = 1, k2 = N / 2;
            int x, y, z, i, j, l;

            for (x = 0; x < log2N; x++)
            {
                l = 0;

                for (y = 0; y < k1; y++, l += k0)
                {
                    for (z = 0; z < k2; z++)
                    {
                        i = z + l; j = i + k2;

                        x_even = data[i];
                        x_odd = data[j];

                        data[i] = x_even + x_odd;
                        data[j] = x_even - x_odd;
                    }
                }

                k0 /= 2; k1 *= 2; k2 /= 2;
            }
        }
        #endregion
    }
}
