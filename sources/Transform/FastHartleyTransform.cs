using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Hartley transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
    /// </remarks>
    [Serializable]
    public class FastHartleyTransform : FloatTransformBase, ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FastFourierTransform FFT;
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
        public override Direction Direction
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
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override float[] Forward(float[] A)
        {
            Complex32[] B = Matrice.ToComplex(A);
            B = FFT.Forward(B);

            int length = A.Length, i;
            float[] Hk = new float[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = B[i].Real - B[i].Imag;
            }

            return Hk;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override float[] Backward(float[] B)
        {
            Complex32[] A = Matrice.ToComplex(B);
            A = FFT.Backward(A);

            int length = B.Length, i;
            float[] Hk = new float[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = A[i].Real + A[i].Imag;
            }

            return Hk;
        }
        #endregion
    }
}
